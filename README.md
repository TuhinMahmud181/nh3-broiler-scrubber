# Ammonia recovery from broiler-house exhaust: process and cost model

MATLAB model for the article:

> Mahmud, T., & Sarker, S. (2026). Recovering ammonia from broiler-house exhaust as ammonium sulphate: a model-based design, cost and uncertainty analysis for three climates. *Biosystems Engineering* (submitted).

The model links the daily ammonia emission and ventilation of a broiler house to a spray scrubber or a packed bed with bypass, the acid and product balance, and the annual costs. It runs the design optimisation, the sensitivity and Sobol uncertainty analysis and the scenarios for three locations: a temperate (United States), a hot-arid (Gulf) and a tropical developing-country case.

## Files

| File | Content |
|---|---|
| `model/nh3_model.m` | One function file with local functions: parameters and location data, house model, spray and packed-bed scrubber models, flock simulation (daily or hourly), capital cost and annual economics, model checks, sizing, grid optimisation, one-at-a-time and global sensitivity analysis, scenarios and the journal figures. Every parameter carries its source as a comment. |
| `model/run_nh3_model.m` | Run script. Executes all analyses in the order of the article, prints every table and saves the figures. |

## Requirements

MATLAB R2024a. No toolboxes are needed.

## How to run

1. In MATLAB, set the Current Folder to `model`.
2. Run `run_nh3_model`.

A full run takes about ten minutes on a desktop computer. For a quick run, set `gsa_N` in section 1 of `run_nh3_model.m` to a smaller value, for example 128. Figures are saved in a new folder `figures` next to `model`.

Single tasks can be called directly, for example:

```matlab
p = nh3_model('params');               % base parameters
e = nh3_model('economics', p, 'us');   % one location: 'us', 'gulf' or 'developing'
```

## Sources in the code comments

"Paper N" in the code comments refers to the numbering of the authors' literature database. The corresponding references in the article are:

| Paper | Reference |
|---|---|
| 1 | Hadlocon et al. (2015a) |
| 3 | Hadlocon et al. (2014a) |
| 6 | Van der Heyden et al. (2016) |
| 10 | Moore et al. (2018) |
| 13 | Melse and Ogink (2005) |
| 14 | Hadlocon et al. (2014b) |
| 15 | Hadlocon et al. (2015c) |
| 17 | Moore et al. (2011) |
| 18 | Wheeler et al. (2006) |
| 19 | Gates et al. (2008) |
| 20 | Miles et al. (2014), listed in the Supplementary Material |
| 21 | Billet and Schultes (1999) |
| 22 | Melse et al. (2006) |
| 27 | Moore et al. (2025) |
| 33 | Miles et al. (2012) |

## Licence

MIT licence, see `LICENSE`.

## How to cite

Please cite the article above. Citation details for the code are in `CITATION.cff`.
