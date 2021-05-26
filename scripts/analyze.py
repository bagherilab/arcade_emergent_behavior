import time
import numpy as np
from .utilities import *

# ANALYSIS UTILITY FUNCTIONS ===================================================

def convert(arr, R):
    """Convert hexagonal coordinates to triangular lattice."""
    L = 2*R - 1
    W = 4*R - 2
    offx = R - 1
    offy = 2*R - 2
    xy = [[u + R - 1 - offx, (w - v) + 2*R - 2 - offy] for u, v, w in arr]
    return list(zip(*xy)), offx, offy, L, W

def get_radius(c):
    """Get radius for given coordinates."""
    u, v, w = c
    return int((abs(u) + abs(v) + abs(w))/2.0)

def get_inds(D, seed, time, H, exclude):
    """Get indices for for seed at selected timepoint."""
    return [(i, j, k)
        for k in range(0, 2*H - 1)
        for i, e in enumerate(D['pop'][seed,time,k,:,:])
        for j, p in enumerate(e) if p not in exclude]

def get_rings(R):
    """Gets ring size for each radius."""
    return [1] + [6*i for i in range(1, R)]

def calculate_statistics(data):
    """Calculate statistics for data."""
    return {
        "mean": np.mean(data, axis=0).tolist(),
        "max": np.max(data, axis=0).tolist(),
        "min": np.min(data, axis=0).tolist(),
        "std": np.std(data, axis=0, ddof=1).tolist()
    }

def calculate_nan_statistics(data):
    """Calculate statistic for data excluding nans."""
    D = np.swapaxes(data, 0, 1)
    unnan = [[e for e in d if not np.isnan(e)] for d in D]
    return {
        "mean": [np.mean(d) if d else np.nan for d in unnan],
        "max": [np.max(d) if d else np.nan for d in unnan],
        "min": [np.min(d) if d else np.nan for d in unnan],
        "std": [np.std(d, ddof=1) if d else np.nan for d in unnan]
    }

# GENERAL METRICS ==============================================================

def get_count(inds):
    """Get number of cells."""
    return len(inds)

def get_volume(D, seed, time, inds):
    """Get total volume of cells."""
    volumes = [D['volume'][seed,time,k,i,p] for i, p, k in inds]
    return np.sum(volumes)

def get_diameter(C, inds):
    """Get diameter of colony."""
    if len(inds) == 0:
        return 0

    layers = list(set([k for i, p, k in inds]))
    sort = [np.array([C[i] for i, p, k in inds if k == layer]) for layer in layers]

    deltas = [[np.max(ma - mi + 1, 0) for ma, mi
        in zip(np.max(coords, axis=0), np.min(coords, axis=0))] for coords in sort]
    diams = [np.mean(delta) for delta in deltas]

    return np.max(diams)

def get_cycle(D, seed, time, inds):
    """Get average cell cycle length."""
    cycles = [D['cycle'][seed,time,k,i,p] for i, p, k in inds]
    cycles = list(filter(lambda x : x != -1, cycles))
    return np.mean(cycles) if len(cycles) > 0 else np.nan

def get_symmetry(C, R, inds):
    """Get symmetry of colony."""
    if len(inds) == 0:
        return np.nan

    layers = list(set([k for i, p, k in inds]))
    coord_sorted = [[tuple(C[i]) for i, p, k in inds if k == layer] for layer in layers]
    all_coord_sets = [set(layer) for layer in coord_sorted]
    unique_coord_sets = [set([tuple(sorted(list(np.abs(c)))) for c in coords]) for coords in all_coord_sets]

    symmetries = []

    for unique_coord_set, all_coord_set in zip(unique_coord_sets, all_coord_sets):
        add = 0
        checked = { unique_coord: False for unique_coord in unique_coord_set }
        deltas = []

        for coord in all_coord_set:
            unique = tuple(sorted(list(np.abs(coord))))

            # only check the symmetry if not yet checked
            if not checked[unique]:
                checked[unique] = True

                sym_coords = get_symmetric(coord) # set of symmetric coordinates
                delta_set = sym_coords - all_coord_set # symmetric coordinates not in full set
                deltas.append(len(delta_set)/len(sym_coords))

        numer = np.sum(deltas)
        denom = len(unique_coord_set) + add
        symmetries.append(1 - numer/denom)

    return np.mean(symmetries)

def get_symmetric(coord):
    """Get list of symmetric coordinates."""
    u, v, w = coord
    return {(u, v, w), (-w, -u, -v), (v, w, u), (-u, -v, -w), (w, u, v), (-v, -w, -u)}

def get_pop(D, seed, time, inds, pop):
    """Get cell populations."""
    return len([1 for i, p, k in inds if D['pop'][seed,time,k,i,p] == pop])

def get_type(D, seed, time, inds, typ):
    """Get cell types."""
    return len([1 for i, p, k in inds if D['type'][seed,time,k,i,p] == typ])

def make_outline(C, R, inds):
    """Get outline of colony."""
    coords = [C[i] for i, p, k in inds]
    xy, offx, offy, L, W = convert(coords, R)
    xy = list(zip(xy[0], xy[1]))

    arr = np.zeros((W, L), dtype=np.uint8)

    lines = []

    [arr.itemset((y + offy, x + offx), 1) for x, y in xy]
    [arr.itemset((y + offy + 1, x + offx), 1) for x, y in xy]
    lines = []

    # Draw left and right segments.
    for j, row in enumerate(arr):
        for i, col in enumerate(row):
            if row[i] == 1:
                if row[i - 1] == 0 or i == 0:
                    lines.append(get_left(i, j, offx, offy))
                if i == len(row) - 1 or row[i + 1] == 0:
                    lines.append(get_right(i, j, offx, offy))

    # Draw up and down segments.
    tarr = np.transpose(arr)
    for i, col in enumerate(tarr):
        for j, row in enumerate(col):
            if col[j] == 1:
                if col[j - 1] == 0 or j == 0:
                    lines.append(get_up_down(i, j, offx, offy))
                if j == len(col) - 1 or col[j + 1] == 0:
                    lines.append(get_up_down(i, j, offx, offy))

    return lines

def get_up_down(i, j, offx, offy):
    """Get up/down outline edge."""
    case = 1 if (i + j)%2 == 0 else 0
    x = i
    y = j + (-1 if (i + j)%2 == 0 else 0)
    return [x - offx, y - offy, case]

def get_left(i, j, offx, offy):
    """Get left outline edge."""
    case = 2 if (i + j)%2 == 0 else 3
    x = i
    y = j + (-1 if (i + j)%2 == 0 else 0)
    return [x - offx, y - offy, case]

def get_right(i, j, offx, offy):
    """Get right outline edge."""
    case = 4 if (i + j)%2 == 0 else 5
    x = i
    y = j + (-1 if (i + j)%2 == 0 else 0)
    return [x - offx, y - offy, case]

# TEMPORAL METRICS =============================================================

def get_temporal_counts(T, N, inds):
    """Get cell counts over time."""
    return [[get_count(inds[i][t])
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_temporal_volumes(D, T, N, inds):
    """Get cell volumes over time."""
    return [[get_volume(D, i, t, inds[i][t])
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_temporal_diameters(T, N, C, inds):
    """Get colony diameter over time."""
    return [[get_diameter(C, inds[i][t])
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_temporal_cycles(D, T, N, inds):
    """Get cell cycles over time."""
    return [[get_cycle(D, i, t, inds[i][t])
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_temporal_growths(T, N, C, inds, t0):
    """Get growth metric over time."""
    diams = get_temporal_diameters(T, N, C, inds)
    return [[np.polyfit(T[t0:t], diams[i][t0:t], 1)[0]
        for t in range(t0 + 2, len(T))]
        for i in range(0,N)]

def get_temporal_symmetries(T, N, C, R, inds):
    """Get symmetry metric over time."""
    return [[get_symmetry(C, R, inds[i][t])
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_temporal_activity(D, T, N, inds, types):
    """Get activity metric over time."""
    typs = get_temporal_types(D, T, N, inds, types)

    active_types = [type for i, type in enumerate(typs) if i in [3, 4]]
    inactive_types = [type for i, type in enumerate(typs) if i in [1, 6]]

    active = np.sum(np.array(active_types), axis=0)
    inactive = np.sum(np.array(inactive_types), axis=0)

    total = np.add(active, inactive)
    total[total == 0] = -1 # set total to -1 to avoid divide by zero error
    activity = np.subtract(np.multiply(np.divide(active, total), 2), 1)
    activity[total == -1] = np.nan # set activity to nan for the case of no active or inactive cells

    return activity

def get_temporal_pops(D, T, N, inds, pops):
    """Get cell populations over time."""
    return [[[get_pop(D, i, t, inds[i][t], pop)
        for t in range(0, len(T))]
        for i in range(0, N)]
        for pop in pops]

def get_temporal_types(D, T, N, inds, types):
    """Get cell types over time."""
    return [[[get_type(D, i, t, inds[i][t], typ)
        for t in range(0, len(T))]
        for i in range(0, N)]
        for typ in types]

# SPATIAL METRICS ==============================================================

def get_spatial_counts(C, N, H, inds):
    """Get cell counts array."""
    arr = np.zeros((2*H - 1, len(C)))
    [[np.add.at(arr, (k, i), 1) for i, p, k in ind] for ind in inds]
    return np.divide(arr, N)

def get_spatial_volumes(D, C, N, H, inds):
    """Get cell volumes array."""
    arr = np.zeros((2*H - 1, len(C)))
    [[np.add.at(arr, (k, i), d['volume'][k,i,p]) for i, p, k in ind] for d, ind in zip(D, inds)]
    return np.divide(arr, N)

def get_spatial_types(D, C, N, H, inds, TYPES):
    """Get cell types arrays."""
    arr = np.zeros((2*H - 1, len(C), len(TYPES)))
    [[np.add.at(arr, (k, i, d['type'][k,i,p]), 1) for i, p, k in ind] for d, ind in zip(D, inds)]
    return np.divide(arr, N)

def get_spatial_pops(D, C, N, H, inds, POPS):
    """Get cell populations arrays."""
    arr = np.zeros((2*H - 1, len(C), len(POPS)))
    [[np.add.at(arr, (k, i, d['pop'][k,i,p]), 1) for i, p, k in ind] for d, ind in zip(D, inds)]
    return np.divide(arr, N)

def get_spatial_outlines(C, N, R, H, inds):
    """Get colony outline arrays."""
    _, offx, offy, L, W = convert(C[0:1], R)
    arr = np.zeros((2*H - 1, L, W, 6))
    for ind in inds:
        layers = list(set([k for i, p, k in ind]))
        sort = [[layer, [(i, p, k) for i, p, k in ind if k == layer]] for layer in layers]
        [np.add.at(arr, (layer, line[0] + offx, line[1] + offy, line[2]), 1) for layer, s in sort for line in make_outline(C, R, s)]
    return np.divide(arr, N)

# GENERAL ANALYSIS =============================================================

def analyze_metrics(D, R, H, T, N, C, POPS, TYPES, outfile, code, exclude=[-1], timepoints=[], seeds=[]):
    """Analyze results for metrics across time."""

    offset = timepoints[0] + 2

    start = time.time()
    inds = [[get_inds(D["agents"], j, i, H, exclude)
        for i in range(0, len(T))]
        for j in range(0, N)]
    end = time.time()
    print(end - start)

    counts = get_temporal_counts(T, N, inds)
    _analyze_metrics(counts, T, f"{outfile}{code}", ".METRICS.COUNTS")

    volumes = get_temporal_volumes(D["agents"], T, N, inds)
    volumes = [[int(v) for v in vol] for vol in volumes]
    _analyze_metrics(volumes, T, f"{outfile}{code}", ".METRICS.VOLUMES")

    cycles = get_temporal_cycles(D["agents"], T, N, inds)
    _analyze_metrics_nan(cycles, T, f"{outfile}{code}", ".METRICS.CYCLES")

    diameters = get_temporal_diameters(T, N, C, inds)
    _analyze_metrics(diameters, T, f"{outfile}{code}", ".METRICS.DIAMETERS")

    types = get_temporal_types(D["agents"], T, N, inds, TYPES)
    _analyze_metrics_list(types, T, TYPES, f"{outfile}{code}", ".METRICS.TYPES")

    pops = get_temporal_pops(D["agents"], T, N, inds, POPS)
    _analyze_metrics_list(pops, T, POPS, f"{outfile}{code}", ".METRICS.POPS")

    growth = get_temporal_growths(T, N, C, inds, offset)
    nan_growth = [np.nan] * (offset + 2)
    _analyze_metrics_nan([nan_growth + grow for grow in growth], T, f"{outfile}{code}", ".METRICS.GROWTH")

    symmetry = get_temporal_symmetries(T, N, C, R, inds)
    _analyze_metrics_nan(symmetry, T, f"{outfile}{code}", ".METRICS.SYMMETRY")

    activity = get_temporal_activity(D["agents"], T, N, inds, TYPES)
    _analyze_metrics_nan(activity, T, f"{outfile}{code}", ".METRICS.ACTIVITY")

def analyze_seeds(D, R, H, T, N, C, POPS, TYPES, outfile, code, exclude=[-1], timepoints=[], seeds=[]):
    """Analyze results for metrics grouped by seed."""

    offset = timepoints[0] + 2

    start = time.time()
    inds = [[get_inds(D["agents"], j, i, H, exclude)
        for i in range(0, len(T))]
        for j in range(0, N)]
    end = time.time()
    print(end - start)

    TT = [(T[t], t) for t in timepoints]

    counts = get_temporal_counts(T, N, inds)
    _analyze_seeds(counts, T, f"{outfile}{code}", TT, ".SEEDS.COUNTS")

    volumes = get_temporal_volumes(D["agents"], T, N, inds)
    volumes = [[int(v) for v in vol] for vol in volumes]
    _analyze_seeds(volumes, T, f"{outfile}{code}", TT, ".SEEDS.VOLUMES")

    diameters = get_temporal_diameters(T, N, C, inds)
    _analyze_seeds(diameters, T, f"{outfile}{code}", TT, ".SEEDS.DIAMETERS")

    cycles = get_temporal_cycles(D["agents"], T, N, inds)
    _analyze_seeds(cycles, T, f"{outfile}{code}", TT, ".SEEDS.CYCLES")

    types = get_temporal_types(D["agents"], T, N, inds, TYPES)
    _analyze_seeds_list(types, T, TYPES, f"{outfile}{code}", TT, ".SEEDS.TYPES")

    pops = get_temporal_pops(D["agents"], T, N, inds, POPS)
    _analyze_seeds_list(pops, T, POPS, f"{outfile}{code}", TT, ".SEEDS.POPS")

    growth = get_temporal_growths(T, N, C, inds, offset)
    nan_growth = [np.nan] * (offset + 2)
    _analyze_seeds([nan_growth + grow for grow in growth], T, outfile + code, TT, ".SEEDS.GROWTH")

    symmetry = get_temporal_symmetries(T, N, C, R, inds)
    _analyze_seeds(symmetry, T, f"{outfile}{code}", TT, ".SEEDS.SYMMETRY")

    activity = get_temporal_activity(D["agents"], T, N, inds, TYPES)
    _analyze_seeds(activity, T, f"{outfile}{code}", TT, ".SEEDS.ACTIVITY")

def analyze_locations(D, R, H, T, N, C, POPS, TYPES, outfile, code, exclude=[-1], timepoints=[], seeds=[]):
    """Analyze results for metrics per location."""

    d = np.take(D["agents"], timepoints, axis=1)
    TT = [T[i] for i in timepoints]

    start = time.time()
    inds = [[get_inds(d, j, i, H, exclude)
        for j in range(0, N)]
        for i in range(0, len(TT))]
    end = time.time()
    print(end - start)

    _pops = ",".join(["POP_" + str(p) for p in POPS])
    _types = ",".join(["TYPE_" + str(t) for t in TYPES])

    counts = [get_spatial_counts(C, N, H, ind) for ind in inds]
    volumes = [get_spatial_volumes(d[:,i,:,:,:], C, N, H, ind) for i, ind in enumerate(inds)]
    pops = [get_spatial_pops(d[:,i,:,:,:], C, N, H, ind, POPS) for i, ind in enumerate(inds)]
    types = [get_spatial_types(d[:,i,:,:,:], C, N, H, ind, TYPES) for i, ind in enumerate(inds)]

    xy, offx, offy, L, W = convert(C, R)

    for i, t in enumerate(TT):
        out = []

        for z, cou, vol, pop, typ in zip(range(1 - H, H), counts[i], volumes[i], pops[i], types[i]):
            joined = [xy[0], xy[1], [z] * len(C), cou, vol] + list(zip(*pop)) + list(zip(*typ))
            out = out + [e for e in zip(*joined) if e[3] != 0]

        header = "x,y,z,COUNT,VOLUME," + _pops + "," + _types + "\n"
        save_csv(f"{outfile}{code}", header, list(zip(*out)), f".LOCATIONS.{format_time(t)}")

def analyze_distribution(D, R, H, T, N, C, POPS, TYPES, outfile, code, exclude=[-1], timepoints=[], seeds=[]):
    """Analyze results for cell state and population distributions."""

    arrt = np.zeros((N, len(T), 2*H - 1, R, len(TYPES)), dtype=np.int32)
    arrp = np.zeros((N, len(T), 2*H - 1, R, len(POPS)), dtype=np.int32)

    start = time.time()
    [np.add.at(arrt, (i, j, k, get_radius(C[m]), e), 1)
        for i, seed in enumerate(D["agents"]['type'])
        for j, time in enumerate(seed)
        for k, layer in enumerate(time)
        for m, loc in enumerate(layer)
        for p, e in enumerate(loc) if D["agents"]['pop'][i,j,k,m,p] not in exclude]
    end = time.time()
    print(end - start)

    start = time.time()
    [np.add.at(arrp, (i, j, k, get_radius(C[m]), e), 1)
        for i, seed in enumerate(D["agents"]['pop'])
        for j, time in enumerate(seed)
        for k, layer in enumerate(time)
        for m, loc in enumerate(layer)
        for p, e in enumerate(loc) if D["agents"]['pop'][i,j,k,m,p] not in exclude]
    end = time.time()
    print(end - start)

    r = np.array(get_rings(R) if len(C[0]) == 3 else X.get_rect_rings(R))

    arrtf = arrt / r[:,None]
    arrpf = arrp / r[:,None]

    arrtf_avg = np.mean(arrtf, axis= 0)
    arrtf_max = np.max(arrtf, axis= 0)
    arrtf_min = np.min(arrtf, axis= 0)
    arrtf_std = np.std(arrtf, axis= 0, ddof=1)

    arrpf_avg = np.mean(arrpf, axis= 0)
    arrpf_max = np.max(arrpf, axis= 0)
    arrpf_min = np.min(arrpf, axis= 0)
    arrpf_std = np.std(arrpf, axis= 0, ddof=1)

    out = [[T[j], r] + ([k - H + 1] if H > 1 else [])
            + [arrtf_avg[j,k,r,t] for t in TYPES]
            + [arrtf_max[j,k,r,t] for t in TYPES]
            + [arrtf_min[j,k,r,t] for t in TYPES]
            + [arrtf_std[j,k,r,t] for t in TYPES]
            + [arrpf_avg[j,k,r,p] for p in POPS]
            + [arrpf_max[j,k,r,p] for p in POPS]
            + [arrpf_min[j,k,r,p] for p in POPS]
            + [arrpf_std[j,k,r,p] for p in POPS]
        for j in range(0, len(T))
        for k in range(0, 2*H - 1)
        for r in range(0, R)]

    metrics = ['_avg', '_max', '_min', '_std']
    _pops = ",".join(["POP_" + str(p) + m for m in metrics for p in POPS])
    _types = ",".join(["TYPE_" + str(t) + m for m in metrics for t in TYPES])

    header = 'time,radius,' + ("height," if H > 1 else "") + _types + ',' + _pops + '\n'
    save_csv(f"{outfile}{code}", header, zip(*out), ".DISTRIBUTION")

def analyze_outlines(D, R, H, T, N, C, POPS, TYPES, outfile, code, exclude=[-1], timepoints=[], seeds=[]):
    """Analyzes results for colony outline."""

    d = np.take(D["agents"], timepoints, axis=1)
    TT = [T[i] for i in timepoints]

    start = time.time()
    inds = [[get_inds(d, j, i, H, exclude)
        for j in range(0, N)]
        for i in range(0, len(TT))]
    end = time.time()

    lines = [get_spatial_outlines(C, N, R, H, ind) for ind in inds]

    xy, offx, offy, L, W = convert(C, R)

    for i, t in enumerate(TT):
        out = []

        for z, lin in zip(range(1 - H, H), lines[i]):
            out = out + [[i - offx, j - offy, z, k, C]
                for i, A in enumerate(lin)
                for j, B in enumerate(A)
                for k, C in enumerate(B) if C != 0]

        header = "x,y,z,DIRECTION,WEIGHT\n"
        save_csv(f"{outfile}{code}", header, list(zip(*out)), f".OUTLINES.{format_time(t)}")

def analyze_concentrations(tar, timepoints, keys, outfile, code, exclude=[]):
    """Analyzes concentration profiles."""

    concentrations = ["glucose", "oxygen", "tgfa"]
    seeds = 50 if len(exclude) == 1 else 20
    radius = 34
    out = {}
    arr = np.zeros((seeds, len(timepoints), radius, len(concentrations)))

    i = 0
    for member in tar.getmembers():
        seed = int(re.findall(r'_([0-9]{2})\.json', member.name)[0])
        json = load_json(member, tar=tar)

        if seed in exclude:
            continue

        for t, tind in enumerate(timepoints):
            for c, conc in enumerate(concentrations):
                concs = json['timepoints'][tind]['molecules'][conc]

                if len(concs) > 1:
                    mid = int((len(concs) - 1)/2)
                    cc = concs[mid]
                    arr[i,t,:,c] = np.array(cc + [np.NaN] * (radius - len(cc)))
                else:
                    arr[i,t,:,c] = np.array(concs[0])

        i = i + 1

    out['_X'] = [x for x in range(0,radius)]
    out['_T'] = [json['timepoints'][tind]["time"] for tind in timepoints]

    for c, conc in enumerate(concentrations):
        out[conc] = {
            "mean": np.mean(arr[:,:,:,c], axis=0).tolist(),
            "std": np.std(arr[:,:,:,c], axis=0, ddof=1).tolist(),
            "min": np.min(arr[:,:,:,c], axis=0).tolist(),
            "max": np.max(arr[:,:,:,c], axis=0).tolist()
        }

    save_json(f"{outfile}{code}", out, ".CONCENTRATIONS")

# ------------------------------------------------------------------------------

def _analyze_metrics(data, T, filename, extension):
    """Save analysis from metrics analysis."""
    out = calculate_statistics(data)
    out['_X'] = T
    save_json(filename, out, extension)

def _analyze_metrics_list(data, T, lst, filename, extension):
    """Save analysis from metrics analysis for lists."""
    out = []
    [out.append(calculate_statistics(data[i])) for i in lst]
    _out = {
        "data": out,
        "_X": T
    }
    save_json(filename, _out, extension)

def _analyze_metrics_nan(data, T, filename, extension):
    """Save analysis from metrics analysis excluding nans."""
    out = calculate_nan_statistics(data)
    out['_X'] = T
    save_json(filename, out, extension)

def _analyze_seeds(data, T, filename, times, extension):
    """Save analysis from seeds analysis."""
    out = []
    for t, i in times:
        out.append({ "_": [x[i] for x in data], "time": t })
    save_json(filename, out, extension)

def _analyze_seeds_list(data, T, lst, filename, times, extension):
    """Save analysis from seeds analysis for lists."""
    _out = []
    for l in lst:
        out = []
        for t, i in times:
            out.append({ "_": [x[i] for x in data[l]], "time": t })
        _out.append(out)
    save_json(filename, _out, extension)
