# Ammonia recovery from broiler-house exhaust: process and cost model

MATLAB code for the article:

> Mahmud, T., & Sarker, S. Recovering ammonia from broiler-house exhaust as ammonium sulphate: a techno-economic and uncertainty analysis for three climates. Submitted to the *Journal of Environmental Management*.

The model follows one broiler house through each flock. It links the daily ammonia emission and ventilation of the house to an acid scrubber with bypass (a spray scrubber or a packed bed), the acid and ammonium sulphate balance, and the annual costs. It is run for three cases: a temperate (United States), a hot-arid (Gulf) and a tropical developing-country case.

The run script reproduces all results of the article:

- checks of the model against published laboratory, field and house data
- flock simulations for the cold, mild and hot season
- annual costs and fertiliser income
- scrubber sizing and design optimisation by grid search
- one-at-a-time sensitivity and Sobol uncertainty analysis
- scenarios: capital cost, operating strategy, airflow–emission link, climate against prices, day–night ventilation, calibrated emission and no extrapolation of the spray correlation
- the figures of the article

## Files

- `model/run_nh3_model.m`: run script. Runs all analyses in the order of the article, prints the tables and saves the figures.
- `model/nh3_model.m`: all model equations and data, with the source of every parameter in the comments.

## Requirements

MATLAB R2024a or later. No toolboxes are needed.

## How to run

1. In MATLAB, set the Current Folder to `model`.
2. Run `run_nh3_model`.

A full run takes about ten minutes on a desktop computer. For a quicker run, set `gsa_N` in section 1 of `run_nh3_model.m` to a smaller value, for example 128. The figures are saved in a new folder `figures` next to `model`.

Single tasks can also be called from the Command Window:

```matlab
p = nh3_model('params');               % base parameters
e = nh3_model('economics', p, 'us');   % one location: 'us', 'gulf' or 'developing'
```

## Licence

MIT Licence. See `LICENSE`.

## Citation

If you use this code, please cite the article above. Citation details are in `CITATION.cff`.
