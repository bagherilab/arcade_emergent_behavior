import numpy as np
from .analyze import *

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

