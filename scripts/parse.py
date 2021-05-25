import tarfile
import numpy as np
from glob import glob
from .utilities import *

# PARSING UTILITY FUNCTIONS ====================================================

def get_files(path, name):
    """Gets list of files from directory."""
    return glob(f"{path}**/{name}*.tar.xz") + glob(f"{path}**/{name}*.json")

def get_coords(R):
    """Get coordinates for given radius."""
    return [[u,v,w]
        for u in range(-R + 1,R)
        for v in range(-R + 1,R)
        for w in range(-R + 1,R)
        if (u + v + w) == 0]

def get_struct(c):
    """Convert cell features into tuple."""
    if c[-1]:
        return (c[1], c[2], np.round(c[4]), np.round(np.mean(c[-1])))
    else:
        return (c[1], c[2], np.round(c[4]), -1)

def extract_setup_fields(jsn):
    """Extracts simulation setup fields."""
    R = jsn["config"]["size"]["radius"]
    H = jsn["config"]["size"]["height"]
    time = [tp["time"] for tp in jsn["timepoints"]]
    pops = [p[0] for p in jsn["config"]["pops"]]
    types = [i for i in range(0,7)]
    return R, H, time, pops, types

def extract_agents_fields(lst, coords, H, N):
    """Extract cell agent fields."""

    # Create empty structured array.
    container = np.empty((2*H - 1, len(coords), N),
        dtype = {
            'names': ['pop', 'type', 'volume', 'cycle'],
            'formats': [np.int8, np.int8, np.int16, np.int16]
        })

    # Set all values in array to -1.
    container[:] = -1

    # Compile entries
    [container.itemset((coord[-1] + H - 1, coords.index(coord[0:-1]), cell[3]), get_struct(cell))
        for coord, cells in lst for cell in cells]

    return container

# GENERAL PARSING ==============================================================

def parse_simulations(name, data_path, result_path, exclude):
    """Parses simulation files."""

    for file in get_files(data_path, name):
        # Create empty arrays.
        container = {
            "agents": [],
            "environments": {
                "glucose": [],
                "oxygen": [],
                "tgfa": []
            }
        }

        if is_tar(file):
            # Parse .tar.xz file.
            tar = tarfile.open(file, "r:xz")
            
            # Iterate through all members of the tar.
            for i, member in enumerate(tar.getmembers()):
                seed = int(member.name.replace(".json", "").split("_")[-1])
                
                # Skip if seed is in exclude list.
                if seed in exclude:
                    continue
                    
                print(f"   > {member.name}")
                parse_simulation(load_json(member, tar=tar), container)
        else:
            # Parse .json file
            parse_simulation(load_json(file), container)

        # Compile data.
        data = {
            "agents": np.array(container['agents']),
            "environments": { x: np.array(container['environments'][x], dtype=np.float16)
                for x in container["environments"].keys() },
            "setup": container["setup"]
        }

        # Pickle results.
        save_path = file.replace(".tar.xz", ".pkl").replace(".json", ".pkl").replace(data_path, result_path)
        pickle.dump(data, open(save_path, "wb"), protocol=4)

def parse_simulation(jsn, container):
    """Parse simulation instance."""

    # Get simulation setup.
    R, H, time, pops, types = extract_setup_fields(jsn)
    coords = get_coords(R)
    N = 6

    # Parse agents.
    container["agents"].append([extract_agents_fields(tp["cells"], coords, H, N) for tp in jsn["timepoints"]])

    # Parse environments.
    container["environments"]["glucose"].append([tp["molecules"]["glucose"] for tp in jsn["timepoints"]])
    container["environments"]["oxygen"].append([tp["molecules"]["oxygen"] for tp in jsn["timepoints"]])
    container["environments"]["tgfa"].append([tp["molecules"]["tgfa"] for tp in jsn["timepoints"]])

    # Add simulation setup to container.
    if not "setup" in container:
        container["setup"] = {
            "radius": R,
            "height": H,
            "time": time,
            "pops": pops,
            "types": types,
            "coords": coords
        }
