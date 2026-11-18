% RUN_NH3_MODEL  Run file for the poultry NH3 recovery model.
%
% How to use:
%   1. In MATLAB, set the Current Folder to this 'model' folder.
%   2. Press Run (or type  run_nh3_model  in the command window).
%
% All equations are in nh3_model.m. Change settings only in section 1 below.

clear; clc; close all;
set(0, 'DefaultAxesFontSize', 10, 'DefaultTextFontSize', 10);

%% 1. Settings -----------------------------------------------------------
p = nh3_model('params');          % default parameters, with their sources

% House and flock (change any value, for example):
% p.n_birds        = 28000;       % birds per house
% p.flock_days     = 42;          % market age, d
% p.litter         = 'new';       % 'used' or 'new'
% p.emission_model = 'quadratic'; % 'linear' or 'quadratic'
% p.ER_slope       = 0.024;       % g NH3/bird/d per day of age (0.024-0.038)
% p.emission_airflow = true;      % less emission at low airflow (Miles et al. 2012)
% p.Q_shape        = 'linear';    % ventilation shape: 'paper17' (measured, default) or 'linear'

% Spray scrubber:
% p.n_modules      = 30;          % number of spray modules
% p.n_stages       = 3;           % spray stages per module (1-3)
% p.P_nozzle       = 0.62;        % nozzle pressure, MPa (0.21-0.62)
% p.run_days       = 28;          % scrubber on only from day 1 to this day ([] = whole flock)

% Packed bed:
% p.bed_area       = 6.0;         % m2 face area
% p.bed_height     = 0.90;        % m of packing
% p.pack.C_V       = 0.167;       % packing mass-transfer constant (0.167-0.35)
% p.pH             = 3.0;         % recirculating liquid
% p.u_liquid       = 10;          % m3/h per m2

% Costs:
% p.flocks_per_year      = 6;     % flocks per house per year
% p.econ.AS_value_factor = 1.0;   % share of granular AS price obtained for the liquid product
% p.econ.labour_h_week   = 2;     % hours of attention per week
% Prices for each location are in the local function location_data in nh3_model.m.

types     = {'spray', 'packed'};              % scrubber types to simulate
seasons   = {'cold', 'mild', 'hot'};
locations = {'us', 'gulf', 'developing'};     % location cases for the cost model
sens_type = 'spray';                          % scrubber type for the sensitivity analysis
gsa_N     = 1024;                             % sample size of the global analysis. Model runs =
                                              % gsa_N*13 per location (1024 -> a few minutes in total).
                                              % 1024 is the value used in the paper; 256 is faster; 0 skips the analysis.
save_figs = true;                             % save all figures in ../figures: vector PDF (journal),
                                              % EMF (for Word, Windows only), PNG 600 dpi and single panels

%% 2. Checks against published data --------------------------------------
val = nh3_model('validate',        p);   % spray model vs Hadlocon et al. 2015c; Hadlocon et al. 2014a; Hadlocon et al. 2015a
pk  = nh3_model('validate_packed', p);   % packed-bed model vs Van der Heyden et al. 2016; Melse and Ogink 2005
chk = nh3_model('inputcheck',      p);   % house input model vs Moore et al. 2011 (+ airflow link, Miles et al. 2012)
bp  = nh3_model('bypasscheck',     p);   % bypass load fraction vs Melse et al. 2006

%% 3. Simulate one flock in each season, for each scrubber type ----------
res = cell(numel(types), numel(seasons));
for t = 1:numel(types)
    p.scrubber_type = types{t};
    for k = 1:numel(seasons)
        res{t,k} = nh3_model('flock', p, seasons{k});
    end
end

%% 4. Summary tables -----------------------------------------------------
rows = { ...
    'Inlet NH3, flock mean (ppmv)',      'C_in_mean',       '%10.1f'; ...
    'Inlet NH3, minimum (ppmv)',         'C_in_min',        '%10.1f'; ...
    'Inlet NH3, maximum (ppmv)',         'C_in_max',        '%10.1f'; ...
    'House airflow, mean (m3/h)',        'Q_mean_m3h',      '%10.0f'; ...
    'House airflow, maximum (m3/h)',     'Q_max_m3h',       '%10.0f'; ...
    'Days with air bypassing scrubber',  'days_bypass',     '%10d';   ...
    'NH3 emitted (kg per flock)',        'NH3_emitted_kg',  '%10.0f'; ...
    'NH3 captured (kg per flock)',       'NH3_captured_kg', '%10.0f'; ...
    'Overall capture (%)',               'capture_pct',     '%10.1f'; ...
    'Ammonium sulphate (kg per flock)',  'AS_kg',           '%10.0f'; ...
    'H2SO4 used, pure (kg per flock)',   'acid_kg',         '%10.0f'; ...
    'Max pressure drop, packing (Pa)',   'dp_max_Pa',       '%10.0f'; ...
    'Extra fan energy (kWh per flock)',  'fan_kWh',         '%10.0f'; ...
    'Days inside correlation range',     'days_in_range',   '%10d'};

for t = 1:numel(types)
    s1 = res{t,1}.summary;
    if strcmp(types{t}, 'spray')
        desc = sprintf('%d modules x %d stages', p.n_modules, p.n_stages);
    else
        desc = sprintf('%.1f m2 face area x %.2f m packing, pH %.1f', p.bed_area, p.bed_height, p.pH);
    end
    fprintf('\n=== FLOCK SIMULATION, %s scrubber: %s ===\n', upper(types{t}), desc);
    fprintf('%d birds, %d days. Scrubber capacity %.0f m3/h, liquid recirculated %.0f L/min\n', ...
            p.n_birds, p.flock_days, s1.Q_cap_m3h, s1.liquid_L_min);
    fprintf('%-34s %10s %10s %10s\n', 'Quantity', seasons{:});
    for r = 1:size(rows, 1)
        fprintf('%-34s', rows{r,1});
        for k = 1:numel(seasons)
            fprintf(rows{r,3}, res{t,k}.summary.(rows{r,2}));
        end
        fprintf('\n');
    end
end
fprintf(['\nNotes\n' ...
  '  NaN = not modeled for that scrubber type.\n' ...
  '  Spray: "days inside correlation range" counts the days on which inlet NH3 (30-300 ppmv)\n' ...
  '  and air velocity (about 2-4 m/s) are within the range tested in Hadlocon et al. 2015c. On other days\n' ...
  '  the spray efficiency is an extrapolation and should not be trusted.\n' ...
  '  Packed: physics-based, so no fitted range, but the result depends strongly on the packing\n' ...
  '  constant C_V (see CHECK 3).\n']);

%% 5. Economics: one house, one year, three locations ---------------------
eco = cell(numel(types), numel(locations));
for t = 1:numel(types)
    p.scrubber_type = types{t};
    for k = 1:numel(locations)
        eco{t,k} = nh3_model('economics', p, locations{k});
    end
end

erows = { ...
    'NH3 captured (% of emitted)',        @(e) e.capture_pct,          '%12.1f'; ...
    'Ammonium sulphate (t/yr)',           @(e) e.AS_t,                 '%12.1f'; ...
    'Capital cost (USD)',                 @(e) e.capital,              '%12.0f'; ...
    '  Capital charge (USD/yr)',          @(e) e.cost.capital,         '%12.0f'; ...
    '  Maintenance (USD/yr)',             @(e) e.cost.maintenance,     '%12.0f'; ...
    '  Sulphuric acid (USD/yr)',          @(e) e.cost.acid,            '%12.0f'; ...
    '  Electricity (USD/yr)',             @(e) e.cost.electricity,     '%12.0f'; ...
    '  Water (USD/yr)',                   @(e) e.cost.water,           '%12.0f'; ...
    '  Labour (USD/yr)',                 @(e) e.cost.labour,          '%12.0f'; ...
    '  Product handling (USD/yr)',        @(e) e.cost.handling,        '%12.0f'; ...
    'Total cost (USD/yr)',                @(e) e.total_cost,           '%12.0f'; ...
    'Fertiliser income (USD/yr)',         @(e) e.revenue,              '%12.0f'; ...
    'NET COST (USD/yr)',                  @(e) e.net_cost,             '%12.0f'; ...
    'Net cost per kg NH3 removed (USD)',  @(e) e.cost_per_kg_NH3,      '%12.2f'; ...
    'Net cost per bird (US cents)',       @(e) 100*e.cost_per_bird,    '%12.2f'; ...
    'Net cost per kg live weight (cents)',@(e) 100*e.cost_per_kg_live, '%12.2f'; ...
    'N recovered per bird (g)',           @(e) e.N_recovered_g_bird,   '%12.1f'; ...
    'Break-even AS price (USD/t)',        @(e) e.breakeven_AS_price,   '%12.0f'; ...
    'Acid cost per kg N captured (USD)',  @(e) e.acid_per_kgN,         '%12.2f'; ...
    'AS income per kg N captured (USD)',  @(e) e.AS_value_per_kgN,     '%12.2f'; ...
    'AS price assumed (USD/t)',           @(e) e.prices.price_AS,      '%12.0f'};

for t = 1:numel(types)
    fprintf('\n=== ECONOMICS, %s scrubber, capacity %.0f m3/h: one house, %d flocks/yr, costs in %d USD ===\n', ...
            upper(types{t}), eco{t,1}.capacity_m3h, p.flocks_per_year, p.econ.cost_year);
    fprintf('%-36s %12s %12s %12s\n', 'Quantity', 'US', 'Gulf', 'Developing');
    for r = 1:size(erows, 1)
        fprintf('%-36s', erows{r,1});
        for k = 1:numel(locations)
            fprintf(erows{r,3}, erows{r,2}(eco{t,k}));
        end
        fprintf('\n');
    end
end
fprintf(['\nNotes on the economics\n' ...
  '  Net cost = total cost - fertiliser income. A positive net cost means the farm pays.\n' ...
  '  Break-even AS price = price of ammonium sulphate at which the scrubber pays for itself.\n' ...
  '  Capture percentages carry the model limits listed above (spray extrapolated, packed C_V).\n' ...
  '  Prices marked ASSUMPTION in location_data (nh3_model.m) are placeholders.\n' ...
  '  Acid cost vs AS income per kg N tests the statement of Moore et al. 2025 that the acid\n' ...
  '  costs more than the recovered N is worth.\n']);

%% 6. Sensitivity: which inputs move the cost most ------------------------
p.scrubber_type = sens_type;
sen = cell(1, numel(locations));
for k = 1:numel(locations)
    sen{k} = nh3_model('sensitivity', p, locations{k});
end

%% 7. Sizing: cost and capture against scrubber capacity ------------------
siz = cell(numel(types), numel(locations));
for t = 1:numel(types)
    p.scrubber_type = types{t};
    for k = 1:numel(locations)
        siz{t,k} = nh3_model('sizing', p, locations{k});
    end
end

%% 8. Design optimization: all design variables together ------------------
opt = cell(numel(types), numel(locations));
for t = 1:numel(types)
    p.scrubber_type = types{t};
    for k = 1:numel(locations)
        opt{t,k} = nh3_model('optimize', p, locations{k});
    end
end
fprintf(['\nNotes on the optimization\n' ...
  '  Each row is the cheapest design on the grid that meets the goal in the first column.\n' ...
  '  Packed designs above %d Pa added pressure (packing + mist eliminator) are excluded (p.dp_limit).\n' ...
  '  Spray results at low air velocity rely on the Hadlocon et al. 2015c equation outside its tested range.\n'], p.dp_limit);

%% 9. Global analysis: all uncertain inputs varied together --------------
if gsa_N > 0
    p.scrubber_type = sens_type;
    gsa = cell(1, numel(locations));
    for k = 1:numel(locations)
        gsa{k} = nh3_model('global', p, locations{k}, gsa_N);
    end
    fprintf(['\nNotes on the global analysis\n' ...
      '  Run at the design set in section 1 (not the optimized one), %s scrubber.\n' ...
      '  "total" indices are the ranking to use. "first order" values are noisier at this sample\n' ...
      '  size (about +/- 0.1) and can come out slightly negative; raise gsa_N to tighten them.\n'], sens_type);
end

%% 10. Scenarios -------------------
p.scrubber_type = sens_type;

% (a) capital cost up to the European level of Moore et al. 2025
cx = cell(1, numel(locations));
for k = 1:numel(locations)
    cx{k} = nh3_model('capex', p, locations{k});
end

% (b) scrubber switched on only in the first part of the flock (Moore et al. 2025)
st = cell(1, numel(locations));
for k = 1:numel(locations)
    st{k} = nh3_model('strategy', p, locations{k});
end

% (c) emission reduced at low airflow (Miles et al. 2012)
q = p;  q.emission_airflow = true;
fprintf('\n=== EMISSION-AIRFLOW LINK (Miles et al. 2012), %s scrubber: base model vs link on ===\n', upper(sens_type));
fprintf('%-34s %12s %12s %12s\n', 'Quantity', 'US', 'Gulf', 'Developing');
e0 = cell(1,3);  e1 = cell(1,3);
for k = 1:numel(locations)
    e0{k} = nh3_model('economics', p, locations{k});
    e1{k} = nh3_model('economics', q, locations{k});
end
lines = {'NH3 emitted, t/yr (base)',      @(e) e.NH3_emitted_kg/1000, e0; ...
         'NH3 emitted, t/yr (link on)',   @(e) e.NH3_emitted_kg/1000, e1; ...
         'Capture % (base)',              @(e) e.capture_pct,         e0; ...
         'Capture % (link on)',           @(e) e.capture_pct,         e1; ...
         'USD per kg NH3 (base)',         @(e) e.cost_per_kg_NH3,     e0; ...
         'USD per kg NH3 (link on)',      @(e) e.cost_per_kg_NH3,     e1; ...
         'cents per bird (base)',         @(e) 100*e.cost_per_bird,   e0; ...
         'cents per bird (link on)',      @(e) 100*e.cost_per_bird,   e1};
for r = 1:size(lines,1)
    fprintf('%-34s', lines{r,1});
    for k = 1:numel(locations), fprintf('%12.2f', lines{r,2}(lines{r,3}{k})); end
    fprintf('\n');
end

% (e) no extrapolation of the spray correlation below 30 ppmv and 2 m/s
q = p;  q.no_extrapolation = true;
fprintf('\n=== NO EXTRAPOLATION OF THE SPRAY CORRELATION, %s scrubber ===\n', upper(sens_type));
fprintf('%-34s %12s %12s %12s\n', 'Quantity', 'US', 'Gulf', 'Developing');
x0 = cell(1,3);  x1 = cell(1,3);
for k = 1:numel(locations)
    x0{k} = nh3_model('economics', p, locations{k});
    x1{k} = nh3_model('economics', q, locations{k});
end
lines = {'Capture % (base)',              @(e) e.capture_pct,         x0; ...
         'Capture % (no extrapolation)',  @(e) e.capture_pct,         x1; ...
         'Income/cost % (base)',          @(e) 100*e.revenue/e.total_cost, x0; ...
         'Income/cost % (no extrap.)',    @(e) 100*e.revenue/e.total_cost, x1; ...
         'USD per kg NH3 (base)',         @(e) e.cost_per_kg_NH3,     x0; ...
         'USD per kg NH3 (no extrap.)',   @(e) e.cost_per_kg_NH3,     x1; ...
         'cents per bird (base)',         @(e) 100*e.cost_per_bird,   x0; ...
         'cents per bird (no extrap.)',   @(e) 100*e.cost_per_bird,   x1};
for r = 1:size(lines,1)
    fprintf('%-34s', lines{r,1});
    for k = 1:numel(locations), fprintf('%12.2f', lines{r,2}(lines{r,3}{k})); end
    fprintf('\n');
end

% (d) further scenarios: climate vs economics, best and worst
%     corner, cost per spray stage, day-night airflow, calibrated emission
rv = nh3_model('further', p);

%% 11. Figures for the manuscript -------------------------------------------
% Every figure is drawn at its printed size (cm), with 8 pt Arial text, so nothing
% is scaled later: what you see on screen is what appears in the journal.
% The drawing code is in nh3_model.m (task 'figures'). Saved in ../figures as
%   vector PDF (submit to the journal), EMF (insert in Word; Windows only),
%   PNG 600 dpi, and each panel on its own (folder ../figures/panels).
% Close any program that has these files open (image viewer, PDF reader) first.
D = struct('res', {res}, 'val', val, 'pk', pk, 'eco', {eco}, 'sen', {sen}, ...
           'siz', {siz}, 'opt', {opt}, 'cx', {cx}, 'st', {st}, ...
           'types', {types}, 'seasons', {seasons}, 'locations', {locations});
if gsa_N > 0, D.gsa = gsa; end
nh3_model('figures', D, save_figs, fullfile(pwd, '..', 'figures'));
