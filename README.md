Supporting code for the article:

> J Yu and N Bagheri. (2020). Agent-based models predict emergent behavior of heterogeneous cell populations in dynamic microenvironments. _Frontiers in Bioengineering and Biotechnology_. doi: [10.3389/fbioe.2020.00249](https://doi.org/10.3389/fbioe.2020.00249)

## Setup files

The `setups` directory contains all the setup files used for running simulations.
Simulations were run using **[ARCADE v2.2](https://github.com/bagherilab/ARCADE/releases/tag/v2.2)**.

The simulations `DEFAULT_random`, `MODULE_COMPLEXITY_both`, `MODULE_COMPLEXITY_metabolism`, and `MODULE_COMPLEXITY_signaling` use modified code; see supplementary materials for details.

## Simulation data

Raw simulation data and results are available on Mendeley Data:

- `DEFAULT` . [http://dx.doi.org/10.17632/wzs7pxkgb9](http://dx.doi.org/10.17632/wzs7pxkgb9)
- `MODULE COMPLEXITY` . [http://dx.doi.org/10.17632/7w2cdsrt87](http://dx.doi.org/10.17632/7w2cdsrt87)
- `PARAMETER SENSITIVITY` . [http://dx.doi.org/10.17632/wpfycmb5sv](http://dx.doi.org/10.17632/wpfycmb5sv)
- `GROWTH CONTEXT` . [http://dx.doi.org/10.17632/mmnh9hsv5y](http://dx.doi.org/10.17632/mmnh9hsv5y)
- `CELL COMPETITION` (1/3) . [http://dx.doi.org/10.17632/5h6ng2y4fc](http://dx.doi.org/10.17632/5h6ng2y4fc)
- `CELL COMPETITION` (2/3) . [http://dx.doi.org/10.17632/gzfwhgdtwz](http://dx.doi.org/10.17632/gzfwhgdtwz)
- `CELL COMPETITION` (3/3) . [http://dx.doi.org/10.17632/dbj46mw6j9](http://dx.doi.org/10.17632/dbj46mw6j9)
- `POPULATION HETEROGENEITY` (1/5) . [http://dx.doi.org/10.17632/bk9s769t3t](http://dx.doi.org/10.17632/bk9s769t3t)
- `POPULATION HETEROGENEITY` (2/5) . [http://dx.doi.org/10.17632/mjpcm5my8z](http://dx.doi.org/10.17632/mjpcm5my8z)
- `POPULATION HETEROGENEITY` (3/5) . [http://dx.doi.org/10.17632/hskf9mxbpw](http://dx.doi.org/10.17632/hskf9mxbpw)
- `POPULATION HETEROGENEITY` (4/5) . [http://dx.doi.org/10.17632/ynzsfdswjz](http://dx.doi.org/10.17632/ynzsfdswjz)
- `POPULATION HETEROGENEITY` (5/5) . [http://dx.doi.org/10.17632/p46dwfyvzg](http://dx.doi.org/10.17632/p46dwfyvzg)

## Pipeline notebooks

#### Parse simulation outputs

The **[`parse_simulation_outputs`](parse_simulation_outputs.ipynb)** notebook provides the functions and scripts for parsing simulation files (`.json`) into pickled numpy arrays (`.pkl`).
These parsed results are included with the raw simulation data.

#### Analyze data & results

The **[`analyze_data_results`](analyze_data_results.ipynb)** notebook provides functions and scripts for running basic analysis on simulation data and parsed results.
All resulting `.json` and `.csv` files are provided in the `analysis` directory.

#### Generate figure inputs

The **[`generate_figure_inputs`](generate_figure_inputs.ipynb)** notebook walks through all the steps necessary to generate figure input files from raw data, parsed files, and basic analysis files.
All resulting files are provided in the `analysis` directory.
Refer to figure section in notebook for more details.

To view figures, start a local HTTP server from the root folder, which can be done using Python or PHP:

```bash
$ python3 -m http.server
$ php -S 127.0.0.1:8000
```

Note that the links in the notebook to figures assume the local port 8000; if your server is running on a different port, the links to the figures from the notebook will not work.
Instead, you can navigate to `http://localhost:XXXX/` where `XXXX` is the port number and follow links to the figures.
