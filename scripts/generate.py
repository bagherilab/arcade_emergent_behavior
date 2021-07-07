import numpy as np
from .analyze import *
from math import sqrt
from scipy.optimize import curve_fit
from scipy.stats import sem

# GENERAL ANALYSIS FUNCTIONS ===================================================

def merge_metrics(file, out, keys, extension, code, tar=None):
    """Merge metrics across conditions."""
    filepath = f"{file}{code}{extension}.json"

    if tar:
        D = load_json(filepath.split("/")[-1], tar=tar)
    else:
        D = load_json(filepath)

    metric = extension.split(".")[-1]
    keys["_Y"] = D['mean']

    keys.pop('time', None)
    out['data'].append(keys)
    out['_X'] = D['_X']

def merge_seeds(file, out, keys, extension, code, tar=None):
    """Merge seed files across conditions."""
    filepath = f"{file}{code}{extension}.json"

    if tar:
        D = load_json(filepath.split("/")[-1], tar=tar)
    else:
        D = load_json(filepath)

    keys.pop('time', None)
    metric = extension.split(".")[-1]

    if metric == "POPS" or metric == "TYPES":
        timepoints = []

        for i, d in enumerate(D):
            for key in keys.keys():
                for dd in d:
                    dd[key] = keys[key]

            if i == 0:
                for dd in d:
                    dd['_'] = [dd['_']]
                timepoints = timepoints + d
            else:
                for tp, dd in zip(timepoints, d):
                    tp['_'].append(dd['_'])
                pass

        out['data'] = out['data'] + timepoints
    else:
        for key in keys.keys():
            for d in D:
                d[key] = keys[key]
        out['data'] = out['data'] + D

def merge_concentrations(file, out, keys, extension, code, tar=None):
    """Merge concentrations across conditions."""
    filepath = f"{file}{code}{extension}.json"

    if tar:
        D = load_json(filepath.split("/")[-1], tar=tar)
    else:
        D = load_json(filepath)

    index = D['_T'].index(unformat_time(keys['time']))
    keys['glucose'] = D['glucose']['mean'][index]
    keys['oxygen'] = D['oxygen']['mean'][index]
    keys['tgfa'] = D['tgfa']['mean'][index]

    out['data'].append(keys)
    out['_X'] = D['_X']

def merge_fractions(file, out, keys, extension, code, tar=None):
    """Merge cell population fractions across conditions."""
    tar_counts = load_tar(file, ".METRICS.COUNTS")
    tar_pops = load_tar(file, ".METRICS.POPS")

    filepath_counts = f"{file}{code}.METRICS.COUNTS.json"
    filepath_pops = f"{file}{code}.METRICS.POPS.json"

    if tar_counts:
        D0 = load_json(filepath_counts.split("/")[-1], tar=tar_counts)
    else:
        D0 = load_json(filepath_counts)

    if tar_pops:
        D = load_json(filepath_pops.split("/")[-1], tar=tar_pops)
    else:
        D = load_json(filepath_pops)

    d0 = D0['mean']
    d = D['data']

    n = len(D['data'])
    F = [np.divide(d[pop]['mean'], d0).tolist() for pop in range(n)]

    keys.pop('time', None)
    for pop in range(n):
        out_pop = keys.copy()
        out_pop['_Y'] = F[pop]
        out_pop["pop"] = pop
        out['data'].append(out_pop)

    out['_X'] = D0['_X']

def merge_maps(file, out, keys, extension, code, tar=None):
    """Merge metrics for heatmap."""
    filepath = f"{file}{code}{extension}.json"

    if tar:
        D = load_json(filepath.split("/")[-1], tar=tar)
    else:
        D = load_json(filepath)

    index = D["_X"].index(keys["time"])
    data = list(keys.values()) + [D["mean"][index]]
    out['data'].append(data)
    out['header'] = list(keys.keys()) + ["metric"]

# ------------------------------------------------------------------------------

def save_metrics(file, extension, out):
    """Save merged metrics files."""
    save_json(file, out, extension)

def save_seeds(file, extension, out):
    """Save merged seed files."""
    save_json(file, out, extension)

def save_concentrations(file, extension, out):
    """Save merged concentrations."""
    save_json(file, out, extension)

def save_fractions(file, extension, out):
    """Save merged fractions."""
    save_json(file, out, extension)

def save_maps(file, extension, out):
    """Save merged metrics map."""
    header = ",".join(out["header"]) + "\n"
    save_csv(file, header, list(zip(*out["data"])), extension.replace("METRICS", "MAPS"))

# DEFAULT STATES ===============================================================

def adjust_portion_size(arr, n):
    """Update portion distribution given max number of portions."""
    if arr:
        while np.sum(arr) > n:
            arr[arr.index(max(arr))] -= 1
        while np.sum(arr) < n:
            arr[arr.index(min(arr))] += 1
    return arr

def get_state_fractions(data, H, C, inds):
    """Get fractions of cells in each cell state."""
    n = 6 # number of triangles in a hex
    arr = np.zeros((2*H - 1, len(C), n))
    arrp = np.empty((2*H - 1, len(C), n))
    arrt = np.empty((2*H - 1, len(C), n))
    arrp[:] = np.nan
    arrt[:] = np.nan

    [np.add.at(arr, (k,i,p), data['volume'][k,i,p])for i, p, k in inds]
    [arrp.itemset((k,i,p), data['pop'][k,i,p]) for i, p, k in inds]
    [arrt.itemset((k,i,p), data['type'][k,i,p]) for i, p, k in inds]

    totals = np.sum(arr, axis=2)
    fracs = [[[j/totals[k][i]*n for j in row if j > 0] for i, row in enumerate(layer)] for k, layer in enumerate(arr)]

    portions = [[[int(round(v)) for v in f] for f in frac] for frac in fracs] # round to integer
    portions = [[[max(1, v) for v in p] for p in portion] for portion in portions] # minimum size
    portions = [[adjust_portion_size(p, n) for p in portion] for portion in portions] # sum to total
    points = [[[0] + np.cumsum(p).tolist()[:-1] if p else [] for p in portion] for portion in portions] # get cumulative counts

    pops = [[[int(j) for j in row if ~np.isnan(j)] for i, row in enumerate(layer)] for k, layer in enumerate(arrp)]
    types = [[[int(j) for j in row if ~np.isnan(j)] for i, row in enumerate(layer)] for k, layer in enumerate(arrt)]

    return portions, points, types, pops

def make_default_states(D, R, H, T, N, C, POPS, TYPES, outfile, code, exclude=[-1], timepoints=[], seeds=[]):
    """Extract cell states for each position."""
    d = np.take(D["agents"], timepoints, axis=1)
    d = np.take(d, seeds, axis=0)
    TT = [T[i] for i in timepoints]

    inds = [[get_inds(d, j, i, H, exclude)
        for j in range(0, len(seeds))]
        for i in range(0, len(TT))]

    _pops = ",".join(["POP_" + str(p) for p in POPS])
    _types = ",".join(["TYPE_" + str(t) for t in TYPES])

    xy, offx, offy, L, W = convert(C, R)

    for i, t in enumerate(TT):
        for j, n in enumerate(seeds):
            ind = [inds[i][j]]
            portions, points, types, pops = get_state_fractions(d[j,i,:,:,:], H, C, ind[0])

            out = []
            for z, portion, point, typ, pop in zip(range(1 - H, H), portions, points, types, pops):
                out = out + [list(i) + [z, a, b, c, d]
                    for i, poi, por, ty, po in zip(zip(*xy), point, portion, typ, pop) if len(por) > 0
                    for j, a, b, c, d in zip([i] * len(poi), poi, por, ty, po)]


            header = "x,y,z,i,n,TYPES,POPS\n"
            save_csv(f"{outfile}{code}", header, list(zip(*out)), f".POSITIONS.{format_time(t)}")

def merge_default_states(file, out, keys, extension, code, tar=None):
    """Merge cell states files across conditions."""
    filepath = f"{file}{code}{extension}.{format_time(keys['time'])}.csv"

    if tar:
        D = load_csv(filepath.split("/")[-1], tar=tar)
    else:
        D = load_csv(filepath)

    d = [[keys['case'], keys['time']] + e for e in D[1:]]
    out['data'] = out['data'] + d
    out['header'] = ["case", "time"] + D[0]

def save_default_states(file, extension, out):
    """Save merged cell states file."""
    save_csv(file, ','.join(out['header']) + "\n", zip(*out['data']), extension)

# DEFAULT SINGLES ==============================================================

def make_default_singles(file, out):
    """Get single seed tracjectories for default simulations."""
    c_file = f"{file}DEFAULT/DEFAULT_C.pkl"
    D, R, H, T, N, C, _, _ = load(c_file)

    c_inds = [[get_inds(D["agents"], j, i, H, [-1]) for i in range(0, len(T))] for j in range(0, N)]
    c_counts = get_temporal_counts(T, N, c_inds)
    c_diams = get_temporal_diameters(T, N, C, c_inds)

    h_file = f"{file}DEFAULT/DEFAULT_H.pkl"
    D, R, H, T, N, C, _, _ = load(h_file)
    h_inds = [[get_inds(D["agents"], j, i, H, [-1]) for i in range(0, len(T))] for j in range(0, N)]
    h_counts = get_temporal_counts(T, N, h_inds)
    h_diams = get_temporal_diameters(T, N, C, h_inds)

    singles = {
        'T': T,
        'COUNTS_C': c_counts,
        'DIAMETERS_C': c_diams,
        'COUNTS_H': h_counts,
        'DIAMETERS_H': h_diams
    }

    save_json(out, singles, f".SINGLES")

# NCI-60 =======================================================================

def calculate_doubling_time(seed, inds, tp):
    """Calculates doubling time between given timepoints."""
    n0 = get_count(inds[seed][0])
    nf = get_count(inds[seed][1])
    return (tp[1]*24 - tp[0]*24)/((np.log(nf) - np.log(n0))/np.log(2))

def calculate_doubling_times(N, tp, inds):
    """Calculates doubling time for all seeds."""
    return [calculate_doubling_time(i, inds, tp) for i in range(0, N)]

def make_doubling_time(file, out):
    """Extract doubling time from simulations."""
    D, R, H, T, N, C, _, _ = load(f"{file}DEFAULT/DEFAULT_C.pkl")
    tp = [2, 16]
    inds = [[get_inds(D["agents"], j, i, H, [-1]) for i in tp] for j in range(0, N)]
    data = calculate_doubling_times(N, [T[tp[0]], T[tp[1]]], inds)

    print(np.mean(data), np.std(data, ddof=1))
    doubling = { "data": data, "mean": np.mean(data) }
    save_json(out, doubling, f".DOUBLING")

def func_exponential(t, no, r):
    """Exponential fit function."""
    return no*np.exp(np.array(t)*r)

def make_exponential_fit(file, out):
    """Fits simulation data with exponential equation."""
    jsn = load_json(f"{out}.SINGLES.json")
    T = jsn['T']
    diams = [[x*30.0 for x in j] for j in jsn['DIAMETERS_C']]
    times = [x for x in T]

    diams = [diam[2:] for diam in diams] # remove first two timepoints from diameters
    times = times[:-2] # adjust time points by 1 day

    fits = [np.polyfit(times, j, 1) for j in diams]
    slopes = [f[0] for f in fits]
    print(np.mean(slopes), np.std(slopes, ddof=1), sem(slopes))

    counts = [j[2:17] for j in jsn['COUNTS_C']]
    t = jsn['T'][0:15]
    fits = [curve_fit(func_exponential, t, y, [0.1, 0.1]) for y in counts]
    r2s = []

    exponential = []

    for y, fit in zip(counts, fits):
        popt, pcov = fit
        residuals = y - func_exponential(t, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r2 = 1 - (ss_res/ss_tot)
        r2s.append(r2)
        exponential.append(np.log(2)/popt[1]*24) # calculate doubling based on exponential fit

    print(np.mean(r2s), np.std(r2s, ddof=1))
    print(np.mean(exponential), np.std(exponential, ddof=1))

    save_json(out, exponential, f".EXPONENTIAL")

# MOON =========================================================================

def make_moon_fit(file, out):
    """Extracts simulation values for fitting moon equation."""
    D, R, H, T, N, C, _, _ = load(f"{file}DEFAULT/DEFAULT_C.pkl")
    inds = [[get_inds(D["agents"], j, i, H, [-1]) for i in range(0, len(T))] for j in range(0, N)]
    counts = get_temporal_counts(T, N, inds)
    volumes = get_temporal_volumes(D["agents"], T, N, inds)

    colony_diameters = get_temporal_diameters(T, N, C, inds)
    average_volumes = [[v/c if c != 0 else np.nan for c, v in zip(cou, vol)] for cou, vol in zip(counts, volumes)]
    cell_diameters = [[2*sqrt(average_volumes[seed][time]/np.pi/4.35)
        for time in range(0, len(T))]
        for seed in range(0, N)]

    elements = zip(*[[n, e, counts[n][t], colony_diameters[n][t], cell_diameters[n][t]]
        for n in range(0, N) for t, e in enumerate(T)])
    header = 'seed,time,count,colony_diameter,cell_diameter\n'
    save_csv(out, header, elements, ".MOON")

def func_moon(X, a, b, c):
    """Equation relating colony diameter, average cell volume, and cell diameter."""
    D, d = X
    return a*np.power(D,b)/np.power(d,c)

def calculate_moon_fit(file, out):
    """Calculate git of moon equation to simulated data."""
    with open(f"{out}.MOON.csv", 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        data = [row for row in reader]

    data = data[1:]
    data = np.array([[float(x) for x in d] for d in data if float(d[1]) >= 1 and float(d[3])*30 < 160])

    N = data[:,2]
    D = data[:,3]*30
    d = data[:,4]

    p0 = 2.4, 2.378, 2.804
    popt, pcov = curve_fit(func_moon, (D, d), N, p0)

    residuals = N - func_moon((D, d), *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((N - np.mean(N))**2)
    r2 = 1 - (ss_res/ss_tot)

    print(popt, r2)

# RANDOM DISTRIBUTION ==========================================================

def make_random_distribution(file, out):
    """Combine cell state distributions for random simulations."""
    tar = load_tar(out, ".DISTRIBUTION")
    values = []

    for case in ["C", "H", "random_R1", "random_R5"]:
        filepath = f"{out}_{case}.DISTRIBUTION.csv"

        if tar:
            D = load_csv(filepath.split("/")[-1], tar=tar)
        else:
            D = load_csv(filepath)

        header = ["case"] + D[0]
        rows = [[case] + r for r in D[1:]]
        values = values + rows

    save_csv(out, ",".join(header) + "\n", zip(*values), ".DISTRIBUTION")

# MODULE CONCENTRATIONS ========================================================

def get_module_concentrations(tar, seeds, timepoints, radius, concentrations):
    """Extracts concentration for all members of given tar."""
    arr = np.zeros((seeds, len(timepoints), radius, len(concentrations)))
    times = []

    for i, member in enumerate(tar.getmembers()):
        json = load_json(member, tar=tar)

        for t, tind in enumerate(timepoints):
            for c, conc in enumerate(concentrations):
                concs = json['timepoints'][tind]['molecules'][conc]
                arr[i,t,:,c] = np.array(concs[0])

            if i == 0:
                times.append(json["timepoints"][tind]["time"])

    return arr, times

def make_module_concentrations(file, out):
    """Extract concentrations from modified module simulations."""

    concentrations = ["glucose", "oxygen", "tgfa"]
    complexities = ["R", "S", "M", "C"]
    timepoints = [2, 4, 16]
    name = out.split("/")[-1]
    radius = 34
    seeds = 20

    for module, code, ext in [("metabolism", "meta", "M"), ("signaling", "sig", "S"), ("both", "", "B")]:
        full_out = []

        for comp1 in complexities:
            for comp2 in complexities:
                if module == "both":
                    keys = { "meta": comp1, "sig": comp2 }
                    filepath = f"{file}{name}/{name}_{module}_{comp1}_{comp2}.tar.xz"
                elif comp2 == "R":
                    keys = { code: comp1 }
                    filepath = f"{file}{name}/{name}_{module}_{comp1}.tar.xz"
                else:
                    continue

                tar = tarfile.open(filepath)
                arr, times = get_module_concentrations(tar, seeds, timepoints, radius, concentrations)

                for t, (tind, tvalue) in enumerate(zip(timepoints, times)):
                    module_out = keys.copy()
                    module_out["time"] = format_time(tvalue)

                    for c, conc in enumerate(concentrations):
                         module_out[conc] = np.mean(arr[:,t,:,c], axis=0).tolist()

                    full_out.append(module_out)

        save_json(out, full_out, f".CONCENTRATIONS.{ext}")

# COMPETITION POPULATIONS ======================================================

def make_competition_populations(D, R, H, T, N, C, POPS, TYPES, outfile, code, exclude=[-1], timepoints=[], seeds=[]):
    """Extract cell populations for each position."""
    d = np.take(D["agents"], timepoints, axis=1)
    d = np.take(d, seeds, axis=0)
    TT = [T[i] for i in timepoints]

    inds = [[get_inds(d, j, i, H, exclude)
        for j in range(0, len(seeds))]
        for i in range(0, len(TT))]

    _pops = ",".join(["POP_" + str(p) for p in POPS])
    _types = ",".join(["TYPE_" + str(t) for t in TYPES])

    xy, offx, offy, L, W = convert(C, R)

    for i, t in enumerate(TT):
        for j, n in enumerate(seeds):
            ind = [inds[i][j]]
            portions, points, types, pops = get_state_fractions(d[j,i,:,:,:], H, C, ind[0])

            out = []
            for z, portion, point, typ, pop in zip(range(1 - H, H), portions, points, types, pops):
                out = out + [list(i) + [z, a, b, c, d]
                    for i, poi, por, ty, po in zip(zip(*xy), point, portion, typ, pop) if len(por) > 0
                    for j, a, b, c, d in zip([i] * len(poi), poi, por, ty, po)]


            header = "x,y,z,i,n,TYPES,POPS\n"
            save_csv(f"{outfile}{code}", header, list(zip(*out)), f".POSITIONS.{format_time(t)}")

def merge_competition_populations(file, out, keys, extension, code, tar=None):
    """Merge cell populations files across conditions."""
    filepath = f"{file}{code}{extension}.{format_time(keys['time'])}.csv"

    if tar:
        D = load_csv(filepath.split("/")[-1], tar=tar)
    else:
        D = load_csv(filepath)

    d = [[keys['param'], keys['perc'], keys['init'], keys['time']] + e for e in D[1:]]
    out['data'] = out['data'] + d
    out['header'] = ["param", "perc", "init", "time"] + D[0]

def save_competition_populations(file, extension, out):
    """Save merged cell populations file."""
    save_csv(file, ','.join(out['header']) + "\n", zip(*out['data']), extension)
