function varargout = nh3_model(task, varargin)
%NH3_MODEL  Poultry-house ammonia recovery model.
%
%   All model equations live in this one file as local functions.
%   Do not run this file directly. Run  run_nh3_model.m  instead.
%
%   Usage (from run_nh3_model.m or the command window):
%     p   = nh3_model('params');                default parameters (a struct)
%     val = nh3_model('validate', p);           check spray model vs published data
%     pk  = nh3_model('validate_packed', p);    check packed-bed model vs published data
%     chk = nh3_model('inputcheck', p);         check house emission input vs published data
%     out = nh3_model('flock', p, 'cold');      simulate one flock ('cold','mild','hot')
%                                               scrubber type is set by p.scrubber_type
%     eco = nh3_model('economics', p, 'us');    annual costs for a location
%                                               ('us', 'gulf' or 'developing')
%     sen = nh3_model('sensitivity', p, 'us');  one-at-a-time sensitivity of the cost
%     siz = nh3_model('sizing', p, 'us');       cost and capture against scrubber size
%     opt = nh3_model('optimize', p, 'us');     best design by exhaustive grid search
%     gsa = nh3_model('global', p, 'us');       Monte Carlo uncertainty + Sobol indices
%     bp  = nh3_model('bypasscheck', p);        bypass load fraction vs Melse et al. 2006 (broilers)
%     st  = nh3_model('strategy', p, 'us');     run the scrubber only in the first part of the flock
%     cx  = nh3_model('capex', p, 'us');        capital-cost scenarios (up to European level, Moore et al. 2025)
%     rv  = nh3_model('further', p);             further scenarios
%     nh3_model('figures', D, true, outdir);    journal figures from the run-file results
%
%   WHAT THE MODEL CONTAINS
%     1. House model   : NH3 emitted per day and exhaust airflow per day
%                        -> NH3 concentration entering the scrubber.
%     2a. Spray scrubber : published semi-mechanistic model, fitted to lab data
%                        (Hadlocon et al. 2015c, Eq. 20 and 21).
%     2b. Packed bed     : physics-based rate model. Gas-film mass transfer
%                        from the Billet & Schultes (1999) correlations, NH3
%                        back-pressure from pH / Henry's law, pressure drop
%                        and extra fan power.
%     3. Mass balance  : NH3 captured, sulphuric acid used, ammonium
%                        sulphate (AS) produced.
%     4. Economics     : capital cost, annual operating cost, fertilizer
%                        income, cost per kg NH3 removed, per bird and per kg
%                        live weight, for three locations.
%
%   WHAT IT DOES NOT CONTAIN YET (later versions)
%     water and heat balance, dust stage, recirculation tank and pH
%     dynamics, crystallizer.
%
%   Paper numbers refer to the '#' column of
%   data/Literature_Data_Extraction.xlsx.

    if nargin < 1
        error('Run run_nh3_model.m, or call nh3_model(''params'').');
    end

    switch lower(task)
        case 'params'
            varargout{1} = default_params();
        case 'validate'
            p = get_params(varargin);
            varargout{1} = run_validation(p);
        case 'validate_packed'
            p = get_params(varargin);
            varargout{1} = run_validation_packed(p);
        case 'inputcheck'
            p = get_params(varargin);
            varargout{1} = run_input_check(p);
        case 'economics'
            p = get_params(varargin);
            if numel(varargin) >= 2
                loc = varargin{2};
            else
                loc = 'us';
            end
            varargout{1} = annual_economics(p, loc);
        case 'sensitivity'
            p = get_params(varargin);
            if numel(varargin) >= 2, loc = varargin{2}; else, loc = 'us'; end
            varargout{1} = run_sensitivity(p, loc);
        case 'optimize'
            p = get_params(varargin);
            if numel(varargin) >= 2, loc = varargin{2}; else, loc = 'us'; end
            varargout{1} = run_optimization(p, loc);
        case 'global'
            p = get_params(varargin);
            if numel(varargin) >= 2, loc = varargin{2}; else, loc = 'us'; end
            if numel(varargin) >= 3, N = varargin{3}; else, N = 256; end
            varargout{1} = run_global(p, loc, N);
        case 'bypasscheck'
            p = get_params(varargin);
            varargout{1} = run_bypass_check(p);
        case 'strategy'
            p = get_params(varargin);
            if numel(varargin) >= 2, loc = varargin{2}; else, loc = 'us'; end
            varargout{1} = run_strategy(p, loc);
        case 'figures'
            % nh3_model('figures', D, save_flag, outdir): journal figures (see run_nh3_model.m)
            make_figures(varargin{:});
        case 'further'
            p = get_params(varargin);
            varargout{1} = run_further(p);
        case 'capex'
            p = get_params(varargin);
            if numel(varargin) >= 2, loc = varargin{2}; else, loc = 'us'; end
            varargout{1} = run_capex(p, loc);
        case 'sizing'
            p = get_params(varargin);
            if numel(varargin) >= 2, loc = varargin{2}; else, loc = 'us'; end
            varargout{1} = run_sizing(p, loc);
        case 'flock'
            p = get_params(varargin);
            if numel(varargin) >= 2
                season = varargin{2};
            else
                season = 'mild';
            end
            varargout{1} = simulate_flock(p, season);
        otherwise
            error('Unknown task "%s". Use params, validate, validate_packed, inputcheck, bypasscheck, flock, economics, sensitivity, sizing, optimize, global, strategy, capex, further or figures.', task);
    end
end

% =====================================================================
%  PARAMETERS
% =====================================================================
function p = default_params()
% Every number has its source. "ASSUMPTION" marks values that are my
% choice, not taken from a paper; change them in run_nh3_model.m.

    % ---- constants
    p.R       = 8.314462;     % J/mol/K
    p.P_atm   = 101325;       % Pa
    p.M_NH3   = 17.031;       % g/mol
    p.M_H2SO4 = 98.079;       % g/mol
    p.M_AS    = 132.14;       % g/mol, (NH4)2SO4

    % ---- house and flock
    p.n_birds      = 25000;   % birds per house (Moore et al. 2011, Table 1: 100,000-122,800 in 4 houses)
    p.flock_days   = 49;      % d, market age (Moore et al. 2011; Wheeler et al. 2006: 42-63 d)
    p.litter       = 'used';  % 'used' (built-up) or 'new' (Wheeler et al. 2006; Gates et al. 2008)
    p.T_exhaust_C  = 25;      % C, exhaust air temperature used outside the flock simulation
    % Exhaust temperature by season (used in the flock simulation). ASSUMPTION
    % (Miles et al. 2014: 20.7-28.9 C in the house; hot-climate houses exhaust warmer air)
    p.T_season = struct('cold', 20, 'mild', 25, 'hot', 31);
    % Day-night swing of the ventilation rate: Q(hour) = Q(day) * (1 + A cos(2 pi (h - 15)/24)),
    % peak at 15:00. A = 0 gives the daily model; A > 0 runs the flock hour by hour.
    p.diurnal_amp  = 0;

    % ---- emission model
    p.emission_model = 'linear';   % 'linear' (Gates et al. 2008, Eq. 1) or 'quadratic' (Miles et al. 2014)
    p.ER_slope       = 0.031;      % g NH3/bird/d per day of age (Wheeler et al. 2006; Gates et al. 2008; farm range 0.024-0.038)
    p.new_litter_lag = 6;          % d with no emission on new litter (Gates et al. 2008, Eq. 1)
    % Optional link between ventilation and emission (Miles et al. 2012):
    % litter releases less NH3 at low airflow. Lab data: airflow 1/3 and 2/3 of the
    % highest gave 59% and 86% of the NH3. Fitted saturating curve
    %     f(x) = (1 - exp(-k x)) / (1 - exp(-k)),  k = 2.3  (gives 0.60 and 0.87)
    % x = house airflow / mild-season airflow on the same day. ASSUMPTION: the
    % emission slope of Wheeler et al. 2006; Gates et al. 2008 is taken as the mild-season value (x = 1).
    p.emission_airflow = false;    % true switches the link on
    p.ea_k             = 2.3;

    % ---- ventilation model
    % Flock-average airflow per house by season (Moore et al. 2011, text of
    % 'Ammonia Emissions in Poultry Houses'), in m3/d:
    p.Q_mean.cold = 0.59e6;   % flocks 4 and 5 (Jan-May): 0.61 and 0.57 million m3/d
    p.Q_mean.mild = 1.31e6;   % flocks 1 and 3: 1.28 and 1.35 million m3/d
    p.Q_mean.hot  = 2.50e6;   % flock 2 (Sept): 2.50 million m3/d
    % Shape over the flock, Q(d) = Q_mean * shape(d), flock mean of shape = 1.
    % 'paper17': daily ventilation of 4 houses, averaged over 5 flocks, digitized
    %            from Moore et al. 2011, Fig. 2A (days 1-49, million m3/d), then divided by
    %            its mean. Low for 3 weeks, rising weeks 4-5, plateau weeks 5-7.
    %            Days 43-49 come from the longer flocks only.
    % 'linear' : Q_mean * (a + b*d/D), a + b/2 = 1 (versions 1-6, an assumption).
    p.Q_shape   = 'paper17';
    p.Q_fig2A   = [0.15 0.25 0.30 0.33 0.30 0.30 0.38 0.38 0.42 0.38 0.45 0.52 0.62 0.58 ...
                   0.68 0.88 0.90 0.83 0.90 0.95 0.92 0.65 0.82 0.98 0.80 1.02 1.30 1.42 ...
                   1.50 1.55 1.45 1.72 1.75 2.10 2.00 1.72 1.75 1.75 2.08 1.62 1.62 1.60 ...
                   1.70 1.80 2.00 2.05 1.90 1.95 2.25];   % million m3/d, days 1-49
    p.Q_shape_a = 0.3;
    p.Q_shape_b = 1.4;

    % ---- which scrubber the flock simulation uses
    p.scrubber_type = 'spray';     % 'spray' or 'packed'

    % ---- packed bed (counter- or cross-current, acid recirculated over packing)
    % Packing: Euroform PN-110 structured plastic packing, the one used in the
    % pig-house scrubber model of Van der Heyden et al. 2016 (Table 1: a, void, C_V, C_L).
    p.pack.a    = 110;       % m2/m3, specific surface area (Van der Heyden et al. 2016, Table 1)
    p.pack.eps  = 0.936;     % void fraction (Van der Heyden et al. 2016, Table 1)
    p.pack.C_V  = 0.167;     % gas-side mass-transfer constant (Billet and Schultes 1999, Table 2b; Van der Heyden et al. 2016, Table 1).
                             % Most other packings in Billet and Schultes 1999 have 0.3-0.45, so this is the
                             % low end; it is the most uncertain packed-bed parameter.
    p.pack.C_h  = 0.511;     % hold-up constant        (Billet and Schultes 1999, Table 2b, Euroform PN-110)
    p.pack.C_P0 = 0.250;     % pressure-drop constant  (Billet and Schultes 1999, Table 2b, Euroform PN-110)
    p.bed_height   = 0.90;   % m of packing in the flow direction (Van der Heyden et al. 2016 reference case)
    p.bed_area     = 6.0;    % m2 face area. ASSUMPTION (design variable)
    p.v_max_packed = 2.0;    % m/s, highest face velocity allowed (Melse and Ogink 2005: 1.7-2.5 m/s)
    p.u_liquid     = 10;     % m3/h of liquid per m2 of bed (Van der Heyden et al. 2016, Table 1)
    p.pH           = 3.0;    % recirculating liquid (Melse and Ogink 2005: 1.3-5; Van der Heyden et al. 2016: 3.2-4)
    p.TAN          = 1.0;    % mol N/L dissolved ammonium in the liquid. ASSUMPTION
                             % (Van der Heyden et al. 2016 discharges at 30 g N/L = 2.1 mol/L)
    p.T_liquid_C   = 20;     % C. ASSUMPTION
    p.fan_eff      = 0.5;    % fan + motor efficiency for the extra pressure drop. ASSUMPTION

    % ---- spray scrubber (laboratory module of Hadlocon et al. 2015c)
    p.col_diameter = 0.356;   % m, round column (Hadlocon et al. 2015c, section 2.4)
    p.z_contact    = 1.65;    % m, contact height used in Eq. 21 (Hadlocon et al. 2015c)
    p.n_stages     = 3;       % spray stages per module, 1 to 3 (Hadlocon et al. 2015c)
    p.P_nozzle     = 0.62;    % MPa, PJ40 nozzle pressure, 0.21 to 0.62 (Hadlocon et al. 2015c, Table 2)
    p.v_max        = 4.0;     % m/s, highest air velocity tested (Hadlocon et al. 2014a; Hadlocon et al. 2015c)
    p.n_modules    = 30;      % number of parallel modules. ASSUMPTION (design variable)
    p.run_days     = [];      % scrubber runs on days 1..run_days only; [] = whole flock.
                              % Moore et al. 2025 suggests running only the first 3-4 weeks.
    % Air above the scrubber capacity (n_modules * area * v_max) bypasses it.

    % ---- utilities and flock schedule (used by the cost model)
    p.dp_spray       = 20;     % Pa, pressure drop of a spray scrubber (Hadlocon et al. 2015a; Hadlocon et al. 2014a; Hadlocon et al. 2014b: 15-26 Pa)
    p.dp_mist        = 20;     % Pa, mist eliminator, added to BOTH scrubber types. ASSUMPTION (typ. 10-30 Pa)
    p.pump_eff       = 0.40;   % pump + motor efficiency. ASSUMPTION (Hadlocon et al. 2015a field unit implies
                               % only about 0.13: 89.5 kWh/d for 49 L/min at 0.59 MPa)
    p.pump_head_packed = 0.15; % MPa, pump pressure for the packed bed (lift + distributor). ASSUMPTION
    p.water_loss_spray  = 1.0; % mL water evaporated per m3 air treated (Hadlocon et al. 2015a: 0.55-1.14)
    p.water_loss_packed = 2.0; % mL/m3 (Van der Heyden et al. 2016: 628 m3/yr at 35,000 m3/h)
    p.flocks_per_year = 6;     % ASSUMPTION: 49 d flock + about 12 d empty (Moore et al. 2011 had 5 flocks/yr)
    p.live_weight_kg  = 2.5;   % kg per bird marketed. ASSUMPTION (Moore et al. 2011: 1.91-3.04 kg)

    % ---- cost basis
    % Capital cost is scaled from published installed costs with
    %   C = C_ref * (capacity / capacity_ref)^n * CEPCI_now / CEPCI_ref
    % Price basis: 2026 USD. Operating prices are 2026 quotes; capital is escalated with the
    % latest annual CEPCI (2025). 2023 = 797.9; 2024 was slightly lower ("just below 800");
    % 2025 was 1.6% above 2024 (Chemical Engineering, 2026) -> about 808. VERIFY in the database.
    p.econ.CEPCI_now  = 808;
    p.econ.cost_year  = 2026;
    % Capital cost of extra spray stages and of packing height. The spray reference unit
    % (Hadlocon et al. 2015a) has three nozzles (stages) per module; each stage fewer saves stage_cost of
    % the capital, each stage more adds it. ASSUMPTION (15%; 0-25% tested in 'further').
    p.econ.stage_cost   = 0.15;
    p.econ.stage_ref    = 3;
    % Packed bed: share of capital that scales with packing height (H_ref = 0.9 m). ASSUMPTION
    p.econ.height_share = 0.30;
    p.econ.height_ref   = 0.9;
    % Liquid product: storage, haulage and spreading of the ammonium sulfate solution.
    % ASSUMPTION: product drawn off at 3 mol N/L (about 200 g AS/L); 5 USD/m3 (range 0-15).
    p.econ.product_N_molL  = 3;
    p.econ.handling_USD_m3 = 5;
    % Spray: Hadlocon et al. 2015a, 14,776 USD for one unit treating 20,727 m3/h (tests 2011-12).
    p.econ.spray_Cref  = 14776;  % USD
    p.econ.spray_Qref  = 20727;  % m3/h
    p.econ.spray_CEPCI = 584.6;  % 2012 annual. VERIFY the index value
    % Exponent from the two published spray units (Hadlocon et al. 2015a: 14,776 USD at 20,727 m3/h;
    % Hadlocon et al. 2014b: 2,645 USD at 1,028 m3/h):  n = ln(14776/2645)/ln(20727/1028) = 0.57
    p.econ.spray_n     = 0.57;
    % Packed: Melse and Ogink 2005, Table 5: 1.3 USD per broiler place for a 30,000-bird house with
    % 75,000 m3/h maximum airflow = 39,000 USD. LOW CONFIDENCE: the paper scaled this
    % from pig houses. Moore et al. 2018 quotes up to 250,000 USD per poultry house in Europe.
    p.econ.packed_Cref  = 39000; % USD
    p.econ.packed_Qref  = 75000; % m3/h
    p.econ.packed_CEPCI = 444.2; % 2004 annual. VERIFY the index value
    p.econ.packed_n     = 0.6;   % six-tenths rule. ASSUMPTION
    p.econ.life_years   = 10;    % Melse and Ogink 2005 (10% depreciation)
    p.econ.maint_frac   = 0.03;  % of capital per year (Melse and Ogink 2005, Table 5)
    p.econ.labour_h_week = 2;    % h/week of attention while birds are in the house. ASSUMPTION
    p.econ.AS_value_factor = 1.0;% fraction of the granular AS price obtained for the liquid
                                 % product. 1.0 follows Hadlocon et al. 2015a; a solution will probably fetch less.

    % ---- multipliers used by the sensitivity analysis (1 = base case)
    p.kya_factor = 1.0;        % on the spray mass-transfer coefficient. The data of Hadlocon et al. 2015c are
                               % 1.17-1.36 times its correlation (CHECK 1), so range 1.0-1.36
    p.Q_scale    = 1.0;        % on the house ventilation rate
    p.econ.capex_factor = 1.0; % on the capital cost
    p.override = struct();     % location prices to replace, e.g. p.override.price_AS = 500
    p.dp_limit = 80;           % Pa (60 packing + 20 mist eliminator), highest TOTAL added pressure (scrubber + mist eliminator) poultry fans tolerate
                               % (Hadlocon et al. 2014a: poultry fans work against 0-60 Pa)
    p.capex_high = 3.0;        % top of the capital-cost factor range in the uncertainty analysis.
                               % Moore et al. 2025: European scrubbers cost over 250,000 USD per house,
                               % about 8x our spray estimate; run_capex shows that case separately.

    % ---- validity range of the Kya correlation (Hadlocon et al. 2015c, below Eq. 21)
    p.range.QL  = [15.31 79.49];    % cm3/s per module
    p.range.D32 = [90 122];         % um
    p.range.C   = [30 300];         % ppmv
    p.range.G   = [81.75 163.50];   % mol/s/m2  (about 2 to 4 m/s)
    p.no_extrapolation = false;     % true: efficiency limited to the edge of the fitted range
end

function p = get_params(args)
    if isempty(args) || isempty(args{1})
        p = default_params();
    else
        p = args{1};
    end
end

% =====================================================================
%  HOUSE MODEL
% =====================================================================
function E = emission_rate(day, p)
% NH3 emitted by the whole house, kg NH3/d, on each day of the flock.
%
% 'linear'    Gates et al. 2008, Eq. 1 (slope measured in Wheeler et al. 2006):
%                ER = 0.031 * x   [g NH3/bird/d]
%                used litter: x = bird age
%                new litter : x = 0 for age < 7, x = age - 6 afterwards
% 'quadratic' Miles et al. 2014: NH3 = 0.0094 d^2 + 0.2542 d + 3.5663  [kg/d per house]
%                measured in a house with about 28,000 birds, so it is
%                scaled here by n_birds/28000.
    day = day(:);
    switch lower(p.emission_model)
        case 'linear'
            if strcmpi(p.litter, 'new')
                x = max(day - p.new_litter_lag, 0);
            else
                x = day;
            end
            E = p.ER_slope .* x .* p.n_birds / 1000;          % g -> kg
        case 'quadratic'
            E = (0.0094*day.^2 + 0.2542*day + 3.5663) * p.n_birds / 28000;
        otherwise
            error('emission_model must be ''linear'' or ''quadratic''.');
    end
end

function Q = ventilation_rate(day, season, p)
% Exhaust airflow of the house, m3/d, on each day of the flock.
    day = day(:);
    season = lower(season);
    if ~isfield(p.Q_mean, season)
        error('season must be ''cold'', ''mild'' or ''hot''.');
    end
    switch lower(p.Q_shape)
        case 'paper17'
            s = p.Q_fig2A(:);
            if p.flock_days ~= numel(s)          % other flock length: stretch the curve
                s = interp1((1:numel(s))'/numel(s), s, day/p.flock_days, 'linear', 'extrap');
                s = max(s, 0.05);
                shape = s / mean(s);
                shape = shape(:);
            else
                shape = s(day) / mean(s);
            end
        case 'linear'
            shape = p.Q_shape_a + p.Q_shape_b * day / p.flock_days;
        otherwise
            error('Q_shape must be ''paper17'' or ''linear''.');
    end
    Q = p.Q_scale * p.Q_mean.(season) * shape;
end

function f = airflow_factor(day, season, p)
% Multiplier on emission from the ventilation-emission link (Miles et al. 2012).
% Returns ones when p.emission_airflow is false.
    day = day(:);
    if ~p.emission_airflow
        f = ones(size(day));  return;
    end
    x = ventilation_rate(day, season, p) ./ ventilation_rate(day, 'mild', setfield(p, 'Q_scale', 1));
    f = (1 - exp(-p.ea_k * x)) / (1 - exp(-p.ea_k));
end

function Vm = molar_volume(T_C, p)
% Ideal-gas molar volume, m3/mol.
    Vm = p.R * (T_C + 273.15) / p.P_atm;
end

function C = concentration_ppmv(E, Q, p)
% NH3 concentration in the exhaust, ppmv, from a mass balance on the house:
%   mol NH3 per day / mol air per day * 1e6
    mol_NH3 = E * 1000 / p.M_NH3;                 % mol/d
    mol_air = Q / molar_volume(p.T_exhaust_C, p); % mol/d
    C = mol_NH3 ./ mol_air * 1e6;
end

% =====================================================================
%  SCRUBBER MODEL  (Hadlocon et al. 2015c)
% =====================================================================
function [QL, D32] = pj40_nozzle(P_MPa)
% One PJ40 full-cone nozzle: liquid flow QL (cm3/s) and Sauter mean droplet
% diameter D32 (um) against operating pressure. Hadlocon et al. 2015c, Table 2
% (manufacturer data), linear interpolation between the three pressures.
    P_tab   = [0.21 0.41 0.62];
    QL_tab  = [15.31 21.64 26.51];
    D32_tab = [122 100 90];
    if P_MPa < P_tab(1) || P_MPa > P_tab(end)
        error('Nozzle pressure must be between 0.21 and 0.62 MPa.');
    end
    QL  = interp1(P_tab, QL_tab,  P_MPa);
    D32 = interp1(P_tab, D32_tab, P_MPa);
end

function Kya = kya_spray(QL, D32, C_in)
% Overall gas-side volumetric mass-transfer coefficient, mol/s/m3 per unit
% mole-fraction difference. Hadlocon et al. 2015c, Eq. 20:
%     Kya = 3916 * QL^0.85 * D32^(-1.30) * C_in^(-0.21)
% QL   total liquid flow in the module, cm3/s
% D32  Sauter mean diameter, um
% C_in inlet NH3, ppmv
    Kya = 3916 .* QL.^0.85 .* D32.^(-1.30) .* C_in.^(-0.21);
end

function eta = scrubber_efficiency(Kya, z, G)
% NH3 removal efficiency (fraction). Hadlocon et al. 2015c, Eq. 17/21:
%     eta = 1 - exp( -Kya * z / G )
% Derivation: plug flow of gas, NH3 reacts instantly and irreversibly with
% the acid, so the driving force is the gas mole fraction itself:
%     G dy/dz = -Kya * y   ->   y_out/y_in = exp(-Kya z / G)
% z  contact height, m
% G  molar gas flux, mol/s/m2
    eta = 1 - exp(-Kya .* z ./ G);
end

function G = gas_molar_flux(v, T_C, p)
% Molar gas flux (mol/s/m2) from superficial air velocity v (m/s).
    G = v ./ molar_volume(T_C, p);
end

function [eta, info] = scrubber_point(C_in, v, P_MPa, n_stages, z, D32_override, QL_override, T_C, p)
% Efficiency of one spray module at one operating point, with range checks.
    [QL1, D32] = pj40_nozzle(P_MPa);
    QL = QL1 * n_stages;                 % each stage has one nozzle (Hadlocon et al. 2015c)
    if ~isempty(D32_override), D32 = D32_override; end
    if ~isempty(QL_override),  QL  = QL_override;  end
    G   = gas_molar_flux(v, T_C, p);
    if p.no_extrapolation
        % Conservative case: no credit from extrapolating the correlation.
        % Concentration and gas flux below the fitted range are raised to the
        % lower edge of the range, so the efficiency never exceeds the value
        % the correlation gives at the edge of its tested conditions.
        C_k = max(C_in, p.range.C(1));
        G_k = max(G, p.range.G(1));
    else
        C_k = C_in;  G_k = G;
    end
    Kya = p.kya_factor * kya_spray(QL, D32, C_k);
    eta = scrubber_efficiency(Kya, z, G_k);
    info.QL = QL;  info.D32 = D32;  info.G = G;  info.Kya = Kya;
    info.in_range = in_range(QL, p.range.QL) & in_range(D32, p.range.D32) & ...
                    in_range(C_in, p.range.C) & in_range(G, p.range.G);
end

function tf = in_range(x, r)
    tol = 1e-3 * max(abs(r));   % 0.1% tolerance for rounding in the published limits
    tf = (x >= r(1) - tol) & (x <= r(2) + tol);
end

% =====================================================================
%  PACKED-BED MODEL  (physics-based)
% =====================================================================
function g = air_properties(T_C, p)
% Dry-air properties and NH3 diffusivity at temperature T_C and 1 atm.
    T = T_C + 273.15;
    g.rho = p.P_atm * 0.028964 / (p.R * T);                 % kg/m3, ideal gas
    g.mu  = 1.716e-5 * (T/273.15)^1.5 * (273.15+110.4)/(T+110.4);   % Pa s, Sutherland
    g.nu  = g.mu / g.rho;                                   % m2/s
    g.D   = 2.28e-5 * (T/298.15)^1.75;                      % m2/s, NH3 in air (0.228 cm2/s at 25 C)
end

function w = water_properties(T_C)
% Liquid properties, taken as water (dilute acid / salt solution).
    w.rho   = 1000 * (1 - (T_C + 288.9414) ./ (508929.2*(T_C + 68.12963)) .* (T_C - 3.9863).^2);  % kg/m3
    w.mu    = 2.414e-5 * 10.^(247.8 ./ (T_C + 273.15 - 140));   % Pa s
    w.nu    = w.mu / w.rho;
    w.sigma = 0.0757 - 1.5e-4 * T_C;                            % N/m
end

function y_eq = nh3_equilibrium_ppmv(pH, TAN, T_C, p)
% NH3 partial pressure above the liquid, in ppmv. This is the "back-pressure":
% the scrubber cannot bring the air below this value.
%   NH4+ <-> NH3(aq) + H+     pKa = 0.09018 + 2729.92/T      (Emerson et al. 1975)
%   NH3(aq) <-> NH3(g)        Henry: kH = 59 mol/L/atm at 25 C, d(ln kH)/d(1/T) = 4200 K
    T      = T_C + 273.15;
    pKa    = 0.09018 + 2729.92 / T;
    f_NH3  = 1 ./ (1 + 10.^(pKa - pH));          % fraction of TAN present as free NH3
    kH     = 59 * exp(4200 * (1/T - 1/298.15));  % mol/L/atm
    p_NH3  = TAN .* f_NH3 ./ kH;                 % atm
    y_eq   = p_NH3 / (p.P_atm/101325) * 1e6;     % ppmv
end

function [kGa, info] = packed_mass_transfer(v, p)
% Volumetric gas-side mass-transfer coefficient kG*a (1/s) of an irrigated
% packing, Billet & Schultes (1999), Trans IChemE 77A, 498-504 (Eq. 6-16).
%
%   hydraulic diameter      d_h = 4*eps/a
%   liquid hold-up          h_L = (12 nu_L u_L a^2 / g)^(1/3) * (a_h/a)^(2/3)
%        a_h/a = C_h Re_L^0.15 Fr_L^0.1            (Re_L < 5)
%        a_h/a = 0.85 C_h Re_L^0.25 Fr_L^0.1       (Re_L >= 5)
%   wetted (interfacial) area
%        a_w/a = 1.5 (a d_h)^-0.5 (u_L d_h/nu_L)^-0.2 (u_L^2 rho_L d_h/sigma)^0.75 (u_L^2/(g d_h))^-0.45
%   gas side
%        kG a_w = C_V (eps - h_L)^-0.5 (a^1.5/d_h^0.5) D (v/(a nu_G))^0.75 (nu_G/D)^(1/3) (a_w/a)
%
% v is the superficial (face) air velocity in m/s.
    g  = 9.81;
    a  = p.pack.a;   eps = p.pack.eps;
    uL = p.u_liquid / 3600;                       % m/s
    G  = air_properties(p.T_exhaust_C, p);
    W  = water_properties(p.T_liquid_C);

    d_h  = 4*eps/a;
    Re_L = uL / (a * W.nu);
    Fr_L = uL^2 * a / g;
    if Re_L < 5
        ah = p.pack.C_h * Re_L^0.15 * Fr_L^0.1;
    else
        ah = 0.85 * p.pack.C_h * Re_L^0.25 * Fr_L^0.1;
    end
    h_L = (12 * W.nu * uL * a^2 / g)^(1/3) * ah^(2/3);

    aw = 1.5 * (a*d_h)^(-0.5) * (uL*d_h/W.nu)^(-0.2) * ...
         (uL^2*W.rho*d_h/W.sigma)^0.75 * (uL^2/(g*d_h))^(-0.45);
    aw = min(aw, 1);                              % cannot wet more than the packing surface

    kGa = p.pack.C_V * (eps - h_L)^(-0.5) * (a^1.5/sqrt(d_h)) * G.D .* ...
          (v/(a*G.nu)).^0.75 * (G.nu/G.D)^(1/3) * aw;

    info.h_L = h_L;  info.aw_over_a = aw;  info.d_h = d_h;  info.Re_L = Re_L;
end

function [eta, info] = packed_point(C_in, v, p)
% Removal efficiency of the packed bed at one operating point.
%
% Gas in plug flow; NH3 that reaches the liquid reacts instantly with the
% acid, so the only resistance is in the gas film, and the liquid exerts a
% constant back-pressure y_eq:
%       v dy/dz = -kG a (y - y_eq)
%   ->  eta = (1 - y_eq/y_in) * (1 - exp(-NTU)),   NTU = kG a H / v
% With y_eq constant the result is the same for counter- and cross-current flow.
    [kGa, info] = packed_mass_transfer(v, p);
    NTU   = kGa .* p.bed_height ./ v;
    y_eq  = nh3_equilibrium_ppmv(p.pH, p.TAN, p.T_liquid_C, p);
    eta   = max(1 - y_eq ./ C_in, 0) .* (1 - exp(-NTU));
    info.kGa = kGa;  info.NTU = NTU;  info.y_eq = y_eq;
end

function dp = packed_pressure_drop(v, p)
% Pressure drop of the irrigated packing, Pa. Billet & Schultes (1999),
% Billet and Schultes 1999, Eq. 28-30, below the loading point, wall effect neglected
% (large bed). The factor (h_L/h_L,S)^0.3 of Eq. 29 is taken as 1, which
% slightly overestimates the pressure drop well below the loading point:
%   dry :  dp0/H = psi0 * a/eps^3 * rho_G v^2 / 2
%          psi0  = C_P0 (64/Re_G + 1.8/Re_G^0.08),  Re_G = v d_p / ((1-eps) nu_G),  d_p = 6(1-eps)/a
%   wet :  dp/dp0 = (eps/(eps - h_L))^1.5 * exp(13300/a^1.5 * sqrt(Fr_L))
% The mist eliminator is NOT included (Hadlocon et al. 2014a measured 15-55 Pa for one).
    a = p.pack.a;  eps = p.pack.eps;
    G = air_properties(p.T_exhaust_C, p);
    [~, info] = packed_mass_transfer(v, p);
    uL   = p.u_liquid / 3600;
    Fr_L = uL^2 * a / 9.81;
    d_p  = 6*(1 - eps)/a;
    Re_G = v .* d_p ./ ((1 - eps) * G.nu);
    psi0 = p.pack.C_P0 * (64 ./ Re_G + 1.8 ./ Re_G.^0.08);
    dp0  = psi0 * a / eps^3 .* G.rho .* v.^2 / 2 * p.bed_height;
    dp   = dp0 * (eps/(eps - info.h_L))^1.5 * exp(13300/a^1.5 * sqrt(Fr_L));
end

% =====================================================================
%  CHECK 3: PACKED-BED MODEL AGAINST PUBLISHED DATA
% =====================================================================
function val = run_validation_packed(p)
% Packed-bed model against one published model result and three full-scale
% measurement programmes. Each case is run twice: with the packing constant
% C_V of Van der Heyden et al. 2016 (0.167) and with a value typical of other packings (0.35),
% because the packings used in the field studies are not identified.

    % {label, a, bed depth (m), face velocity (m/s), pH, inlet ppmv, measured %, note}
    cs = { ...
      'Van der Heyden 2016 (model), pigs',     110, 0.90, 1.62, 4.0, 14.0, 84, 'model vs model'; ...
      'Hol 1999, 30,000 broilers',             150, 0.90, 1.39, 4.0, 18.8, 95, 'depth assumed 0.9 m'; ...
      'Vrielink 1997, 66 pigs',                100, 1.53, 0.93, 3.0,  8.2, 91, 'depth from volume/area'; ...
      'Verdoes 1999, 54 pigs',                 100, 0.90, 0.67, 4.0, 15.6, 99, 'depth from volume/area'};
    CV = [0.167 0.35];

    fprintf('\n=== CHECK 3: packed-bed model (Billet & Schultes) against published data ===\n');
    fprintf('%-40s %9s %12s %12s   %s\n', 'Case', 'Measured', 'C_V = 0.167', 'C_V = 0.35', 'Note');
    n = size(cs,1);
    val.label = cs(:,1);  val.measured = zeros(n,1);  val.model = zeros(n,2);
    for k = 1:n
        q = p;
        q.pack.a     = cs{k,2};
        q.bed_height = cs{k,3};
        q.pH         = cs{k,5};
        q.T_exhaust_C = 22;
        e = zeros(1,2);
        for j = 1:2
            q.pack.C_V = CV(j);
            e(j) = 100 * packed_point(cs{k,6}, cs{k,4}, q);
        end
        val.measured(k) = cs{k,7};  val.model(k,:) = e;
        fprintf('%-40s %9.0f %12.1f %12.1f   %s\n', cs{k,1}, cs{k,7}, e(1), e(2), cs{k,8});
    end
    fprintf(['Removal efficiency in %%. Face velocities are average airflow / face area, with the\n' ...
             'face area from maximum airflow / maximum velocity (Melse and Ogink 2005, Table 2). Liquid load\n' ...
             '10 m3/h/m2 assumed in every case. Void fraction and C_h kept at the Van der Heyden et al. 2016 values.\n']);

    % pressure drop at the reference case
    q = p;  q.T_exhaust_C = 22;
    fprintf('Pressure drop at the Van der Heyden et al. 2016 reference case (0.9 m, 1.62 m/s): %.0f Pa (packing only).\n', ...
            packed_pressure_drop(1.62, q));
    fprintf('NH3 back-pressure of the liquid at 20 C, 1 mol N/L:  pH 3: %.3f   pH 5: %.2f   pH 6: %.1f   pH 7: %.0f ppmv\n', ...
            nh3_equilibrium_ppmv(3, 1, 20, p), nh3_equilibrium_ppmv(5, 1, 20, p), ...
            nh3_equilibrium_ppmv(6, 1, 20, p), nh3_equilibrium_ppmv(7, 1, 20, p));
end

% =====================================================================
%  CHECK 1: SPRAY MODEL AGAINST PUBLISHED DATA
% =====================================================================
function val = run_validation(p)
% Compares model predictions with numbers printed in the papers.
% Each row: description, measured value, model value.

    T = 25;   % C, laboratory air. ASSUMPTION (not stated for each test)
    c = {};   % {label, quantity, measured, predicted, in_range, note}

    % --- A. Mass-transfer coefficient, Hadlocon et al. 2015c section 3.2.4 (30 ppmv, 0.62 MPa)
    [QL1, D32] = pj40_nozzle(0.62);
    c(end+1,:) = {'Hadlocon 2015c: Kya, 1 stage, 30 ppmv, 0.62 MPa', 'Kya', 121.34, ...
                  kya_spray(QL1, D32, 30), true, 'mol/s/m3'};
    c(end+1,:) = {'Hadlocon 2015c: Kya, 3 stages, 30 ppmv, 0.62 MPa', 'Kya', 265.87, ...
                  kya_spray(3*QL1, D32, 30), true, 'mol/s/m3'};

    % --- B. Removal efficiency, Hadlocon et al. 2014a abstract (lab, same type of module)
    [e, i] = scrubber_point(100, 4.0, 0.62, 3, p.z_contact, [], [], T, p);
    c(end+1,:) = {'Hadlocon 2014a: 3 stages, 100 ppmv, 4 m/s, 0.62 MPa', 'eta %', 86, 100*e, i.in_range, ''};
    [e, i] = scrubber_point(400, 4.0, 0.62, 3, p.z_contact, [], [], T, p);
    c(end+1,:) = {'Hadlocon 2014a: 3 stages, 400 ppmv, 4 m/s, 0.62 MPa', 'eta %', 74, 100*e, i.in_range, 'C above 300'};
    [e, i] = scrubber_point(30, 4.0, 0.51, 1, p.z_contact, [], [], T, p);
    c(end+1,:) = {'Hadlocon 2014a: 1 stage, 30 ppmv, 4 m/s, 0.51 MPa', 'eta %', 91, 100*e, i.in_range, ''};
    [e, i] = scrubber_point(30, 3.0, 0.51, 1, p.z_contact, [], [], T, p);
    c(end+1,:) = {'Hadlocon 2014a: 1 stage, 30 ppmv, 3 m/s, 0.51 MPa', 'eta %', 95, 100*e, i.in_range, ''};

    % --- C. Field scrubber, Hadlocon et al. 2015a Table 3 (poultry manure composting).
    % 15 hexagonal modules, 0.16 m2 each, 1.2 m contact height, 3 PJ40 nozzles
    % per module, droplet 113 um. Seasonal inlet and efficiency from Table 3.
    % The airflow during the tests is ambiguous in the paper: 20,727 m3/h
    % with the scrubber, 6,802 m3/h with scrubber and air filter. Both shown.
    A15   = 15 * 0.16;                       % m2
    QLmod = [44.27 49.07 48.20 54.11] * 1000/60 / 15;   % cm3/s per module
    Cin   = [95.15 99.12 92.27 103.63];
    meas  = [75.30 76.78 71.32 80.53];
    seas  = {'winter','spring','summer','autumn'};
    Qair  = [20727 6802];                    % m3/h
    tag   = {'no filter', 'with filter'};
    for k = 1:2
        v = Qair(k) / 3600 / A15;
        for j = 1:4
            [e, i] = scrubber_point(Cin(j), v, 0.59, 3, 1.2, 113, QLmod(j), T, p);
            c(end+1,:) = {sprintf('Hadlocon 2015a: field %s, %s (%.2f m/s)', seas{j}, tag{k}, v), ...
                          'eta %', meas(j), 100*e, i.in_range, ''}; %#ok<AGROW>
        end
    end

    % --- print
    fprintf('\n=== CHECK 1: spray model (Hadlocon et al. 2015c, Eq. 20-21) against published data ===\n');
    fprintf('%-56s %-6s %9s %9s %8s  %s\n', 'Case', 'Qty', 'Measured', 'Model', 'Diff', 'In range');
    n = size(c,1);
    val.label = c(:,1);  val.measured = zeros(n,1);  val.model = zeros(n,1);
    val.in_range = false(n,1);
    for k = 1:n
        val.measured(k) = c{k,3};  val.model(k) = c{k,4};  val.in_range(k) = c{k,5};
        d = c{k,4} - c{k,3};
        if c{k,5}, r = 'yes'; else, r = 'NO (extrapolated)'; end
        fprintf('%-56s %-6s %9.1f %9.1f %+8.1f  %s\n', c{k,1}, c{k,2}, c{k,3}, c{k,4}, d, r);
    end
    fprintf(['Diff = model - measured. For "eta %%" rows it is in percentage points.\n' ...
             'Not checked: Hadlocon et al. 2014b (swine) uses an AAP01 nozzle with 250-329 um droplets,\n' ...
             'outside the 90-122 um range of the correlation.\n']);
end

% =====================================================================
%  CHECK 2: HOUSE INPUT MODEL AGAINST PUBLISHED DATA
% =====================================================================
function chk = run_input_check(p)
% Flock-average NH3 concentration and emission predicted by the house model
% against the measurements of Moore et al. 2011 (four Arkansas broiler houses).
    seasons  = {'cold', 'mild', 'hot'};
    C_meas   = [42.5 20.0 11.1];   % ppm: flocks 4-5 (42.9, 42.2); flocks 1, 3 (19.3, 20.6); flock 2
    E_meas   = [13.6 16.3 14.0];   % kg/d: flocks 4-5 (12.5, 14.7); flocks 1, 3 (17.0, 15.6); flock 2
    day = (1:p.flock_days)';
    fprintf('\n=== CHECK 2: house input model against Moore et al. 2011 (flock averages) ===\n');
    fprintf('%-6s %12s %10s %10s %12s %10s %10s\n', 'Season', 'C meas', 'C model', 'C airflow', 'E meas', 'E model', 'E airflow');
    fprintf('%-6s %12s %10s %10s %12s %10s %10s\n', '', '(ppm)', '', 'link on', '(kg/d)', '', 'link on');
    chk.season = seasons;  chk.C_meas = C_meas;  chk.E_meas = E_meas;
    chk.C_model = zeros(1,3);  chk.E_model = zeros(1,3);  chk.C_link = zeros(1,3);  chk.E_link = zeros(1,3);
    q = p;  q.emission_airflow = true;
    for k = 1:3
        p.T_exhaust_C = p.T_season.(seasons{k});
        Q  = ventilation_rate(day, seasons{k}, p);
        E0 = emission_rate(day, p);
        E1 = E0 .* airflow_factor(day, seasons{k}, q);
        chk.C_model(k) = mean(concentration_ppmv(E0, Q, p));   chk.E_model(k) = mean(E0);
        chk.C_link(k)  = mean(concentration_ppmv(E1, Q, p));   chk.E_link(k)  = mean(E1);
        fprintf('%-6s %12.1f %10.1f %10.1f %12.1f %10.1f %10.1f\n', seasons{k}, C_meas(k), chk.C_model(k), ...
                chk.C_link(k), E_meas(k), chk.E_model(k), chk.E_link(k));
    end
    rmse = @(a,b) sqrt(mean((a-b).^2));
    fprintf('RMSE of emission: %.1f kg/d without the airflow link, %.1f kg/d with it (Miles et al. 2012).\n', ...
            rmse(chk.E_model, E_meas), rmse(chk.E_link, E_meas));
    fprintf(['Measured values are averages of the flocks in that season. The ventilation\n' ...
             'averages are taken from the same paper, so only the emission model and the\n' ...
             'shape over the flock are being tested here. "link on" = emission reduced at low\n' ...
             'airflow (p.emission_airflow = true).\n']);
end

% =====================================================================
%  FLOCK SIMULATION
% =====================================================================
function out = simulate_flock(p, season)
% Simulation of one flock. Quasi-steady: each time step (a day, or an hour
% when p.diurnal_amp > 0) is a steady state at that step's airflow and
% concentration. Exhaust temperature follows the season (p.T_season).

    if isfield(p, 'T_season') && isfield(p.T_season, season)
        p.T_exhaust_C = p.T_season.(season);
    end
    day = (1:p.flock_days)';   D = numel(day);
    E   = emission_rate(day, p) .* airflow_factor(day, season, p);   % kg NH3/d leaving the house
    Q   = ventilation_rate(day, season, p);      % m3/d (daily mean)

    % --- scrubber capacity
    is_packed = strcmpi(p.scrubber_type, 'packed');
    if is_packed
        A_tot = p.bed_area;                          % m2 face area
        v_cap = p.v_max_packed;
    elseif strcmpi(p.scrubber_type, 'spray')
        A_tot = p.n_modules * pi/4 * p.col_diameter^2;
        v_cap = p.v_max;
    else
        error('scrubber_type must be ''spray'' or ''packed''.');
    end
    Q_cap = A_tot * v_cap * 86400;               % m3/d the scrubber can pass
    if isempty(p.run_days), run_days = p.flock_days; else, run_days = min(p.run_days, p.flock_days); end
    on = day <= run_days;                        % scrubber switched on

    % --- time steps: 1 per day, or 24 per day with a day-night swing of the airflow
    A = 0;  if isfield(p, 'diurnal_amp'), A = p.diurnal_amp; end
    if A > 0, nh = 24; else, nh = 1; end
    hr  = repmat((0:nh-1)', D, 1);
    Qx  = kron(Q, ones(nh,1));
    if nh > 1, Qx = Qx .* (1 + A*cos(2*pi*(hr - 15)/24)); end
    Ex  = kron(E, ones(nh,1));                   % emission rate constant within a day
    onx = kron(on, ones(nh,1));
    Cx  = concentration_ppmv(Ex, Qx, p);         % ppmv at scrubber inlet
    Qtx = min(Qx, Q_cap) .* onx;                 % air through the scrubber (rate, m3/d)
    ftx = Qtx ./ Qx;                             % fraction of exhaust treated
    vx  = max(Qtx, 1e-6) / 86400 / A_tot;        % m/s in the scrubber (kept > 0 when off)

    % --- efficiency of the treated stream (all steps at once)
    if is_packed
        etax = packed_point(Cx, vx, p);
        dpx  = packed_pressure_drop(vx, p);      % Pa, packing only
        okx  = onx > 0;                          % physics-based, no fitted range
        dp_fan = dpx + p.dp_mist;                % Pa seen by the fans
    else
        [etax, info] = scrubber_point(Cx, vx, p.P_nozzle, p.n_stages, ...
                                      p.z_contact, [], [], p.T_exhaust_C, p);
        okx  = info.in_range & (onx > 0);
        dpx  = nan(size(vx));                    % not modelled (measured p.dp_spray used)
        dp_fan = (p.dp_spray + p.dp_mist) * ones(size(vx));
    end
    capx = Ex .* ftx .* etax;                    % kg/d captured (rate)
    fanx = Qtx/86400 .* dp_fan / p.fan_eff * 24 / 1000;   % kWh/d (rate)

    % --- back to daily values (means over the hours of a day)
    dm = @(x) mean(reshape(x, nh, D), 1)';
    NH3_cap = dm(capx);
    Q_treat = dm(Qtx);
    f_treat = Q_treat ./ Q;
    C   = dm(Cx);
    v   = dm(vx);
    eta = NH3_cap ./ max(E .* dm(ftx), eps);
    fan_kWh = dm(fanx);
    ok  = all(reshape(okx, nh, D), 1)';
    if is_packed, dp = max(reshape(dpx, nh, D), [], 1)'; else, dp = nan(D,1); end

    % --- mass balance.  2 NH3 + H2SO4 -> (NH4)2SO4
    acid = NH3_cap * p.M_H2SO4 / (2*p.M_NH3);    % kg/d pure H2SO4
    AS   = NH3_cap * p.M_AS    / (2*p.M_NH3);    % kg/d ammonium sulfate
    if is_packed
        L_total = p.u_liquid * p.bed_area * 1000/60;            % L/min recirculated
    else
        [QL1, ~] = pj40_nozzle(p.P_nozzle);
        L_total = QL1 * p.n_stages * p.n_modules * 60/1000;     % L/min recirculated
    end

    out.season  = season;
    out.day     = day;
    out.E       = E;            % kg NH3/d emitted
    out.Q_m3h   = Q / 24;       % m3/h house exhaust
    out.C_in    = C;            % ppmv
    out.f_treat = f_treat;      % fraction of air treated
    out.v       = v;            % m/s in modules
    out.eta     = eta;          % efficiency of treated air
    out.capture = f_treat .* eta;   % overall capture of house NH3
    out.NH3_cap = NH3_cap;      % kg/d
    out.acid    = acid;         % kg/d
    out.AS      = AS;           % kg/d
    out.in_range = ok;          % true where the correlation is inside its tested range
    out.dp      = dp;           % Pa, packing pressure drop (packed bed only)
    out.fan_kWh = fan_kWh;      % kWh/d extra fan energy (scrubber + mist eliminator)

    s.NH3_emitted_kg  = sum(E);
    s.NH3_captured_kg = sum(NH3_cap);
    s.capture_pct     = 100 * sum(NH3_cap) / sum(E);
    s.AS_kg           = sum(AS);
    s.acid_kg         = sum(acid);
    s.C_in_mean       = mean(C);
    s.C_in_min        = min(C);
    s.C_in_max        = max(C);
    s.Q_mean_m3h      = mean(Q)/24;
    s.Q_max_m3h       = max(Q)/24;
    s.Q_cap_m3h       = Q_cap/24;
    s.days_in_range   = sum(ok);
    s.liquid_L_min    = L_total;
    s.dp_max_Pa       = max(dp(on));          % packing only; fans see this + p.dp_mist
    s.days_bypass     = sum(any(reshape(ftx < 1 - 1e-9 & onx > 0, nh, D), 1));
    s.fan_kWh         = sum(fan_kWh);
    s.air_treated_m3  = sum(Q_treat);
    s.run_days        = run_days;
    s.load_treated_pct = 100 * sum(E .* f_treat) / sum(E);   % share of the NH3 load sent through the scrubber
    out.summary = s;
end

% =====================================================================
%  ECONOMICS
% =====================================================================
function L = location_data(name)
% Prices and climate for the three location cases. Prices are 2026 market
% reports found on 15 Nov 2026 (see data/Literature_Data_Extraction.xlsx,
% sheet Product_Economics, for the sources). "ASSUMPTION" = my choice.
%
% season_mix = share of the year's flocks grown under the cold / mild / hot
% ventilation profile of Moore et al. 2011.
    switch lower(name)
        case 'us'
            L.label      = 'United States';
            L.season_mix = [0.4 0.4 0.2];  % Moore et al. 2011: of 5 flocks, 2 cold, 2 mild, 1 hot
            L.price_AS   = 400;    % USD/t ammonium sulphate (2026 reports: 288-474 bulk; 650 retail)
            L.price_acid = 170;    % USD/t H2SO4 (2026 reports: 117-177; one report 318)
            L.price_elec = 0.11;   % USD/kWh. Between US industrial (0.089) and commercial (0.139). ASSUMPTION
            L.price_water = 1.0;   % USD/m3. ASSUMPTION
            L.wage       = 18;     % USD/h farm labour. ASSUMPTION
            L.interest   = 0.06;   % Melse and Ogink 2005 used 6%
            L.range_AS   = [288 650];      % bulk low ... retail
            L.range_acid = [117 318];
            L.range_elec = [0.089 0.139];  % industrial ... commercial average
        case 'gulf'
            L.label      = 'Gulf (Saudi Arabia)';
            L.season_mix = [0 1/3 2/3];    % ASSUMPTION: no cold season. Needs local ventilation data
            L.price_AS   = 250;    % USD/t. ASSUMPTION from Asian export prices (161-278); no Gulf quote found
            L.price_acid = 440;    % USD/t, spot ex-Jeddah, Q3 2026 (unusually high; earlier years far lower)
            L.price_elec = 0.061;  % USD/kWh = 23 halalas, agricultural tariff above 6,000 kWh/month from Jan 2027
            L.price_water = 2.0;   % USD/m3. ASSUMPTION
            L.wage       = 6;      % USD/h. ASSUMPTION
            L.interest   = 0.06;   % ASSUMPTION
            L.range_AS   = [160 330];      % China bulk ... India
            L.range_acid = [170 440];      % North America level ... Jeddah spot Q3 2026
            L.range_elec = [0.048 0.061];  % 18 ... 23 halalas
        case 'developing'
            L.label      = 'Developing country (South Asia)';
            L.season_mix = [0 0.4 0.6];    % ASSUMPTION: tropical, closed commercial houses only
            L.price_AS   = 325;    % USD/t (India, Sept 2026 report; South-East Asia 278)
            L.price_acid = 320;    % USD/t. ASSUMPTION from North-East Asia prices (295-330)
            L.price_elec = 0.105;  % USD/kWh, Bangladesh industrial 11.56-16.06 Tk/kWh (June 2026 tariff)
            L.price_water = 0.3;   % USD/m3. ASSUMPTION
            L.wage       = 1.5;    % USD/h. ASSUMPTION
            L.interest   = 0.10;   % ASSUMPTION (higher cost of borrowing)
            L.range_AS   = [278 400];
            L.range_acid = [170 330];
            L.range_elec = [0.094 0.130];  % off-peak ... peak
        otherwise
            error('location must be ''us'', ''gulf'' or ''developing''.');
    end
end

function [C, Q_cap] = capital_cost(p)
% Installed capital cost in USD of p.econ.cost_year, and capacity in m3/h.
    e = p.econ;
    if strcmpi(p.scrubber_type, 'packed')
        Q_cap = p.bed_area * p.v_max_packed * 3600;
        C = e.capex_factor * e.packed_Cref * (Q_cap/e.packed_Qref)^e.packed_n * e.CEPCI_now/e.packed_CEPCI;
        C = C * ((1 - e.height_share) + e.height_share * p.bed_height / e.height_ref);   % packing height
    else
        Q_cap = p.n_modules * pi/4 * p.col_diameter^2 * p.v_max * 3600;
        C = e.capex_factor * e.spray_Cref * (Q_cap/e.spray_Qref)^e.spray_n * e.CEPCI_now/e.spray_CEPCI;
        C = C * (1 + e.stage_cost * (p.n_stages - e.stage_ref));                         % spray stages
    end
end

function eco = annual_economics(p, loc)
% One year of operation of one house at a location.
%
%   annual quantities = flocks per year * sum over seasons( season share * flock result )
%   capital charge    = capital * CRF,  CRF = i (1+i)^n / ((1+i)^n - 1)
%   net cost          = capital charge + maintenance + acid + electricity
%                       + water + labour - fertilizer income
    L = location_data(loc);
    if isfield(p, 'override') && isstruct(p.override)       % prices replaced by the user
        f = fieldnames(p.override);                          % or by the sensitivity analysis
        for j = 1:numel(f), L.(f{j}) = p.override.(f{j}); end
    end
    seasons = {'cold', 'mild', 'hot'};
    is_packed = strcmpi(p.scrubber_type, 'packed');

    % --- weighted flock results
    q = struct('emit',0,'cap',0,'AS',0,'acid',0,'fan',0,'air',0);
    L_min = 0;  dp_max = 0;
    for k = 1:3
        w = L.season_mix(k);
        if w == 0, continue; end
        o = simulate_flock(p, seasons{k});
        s = o.summary;
        q.emit = q.emit + w * s.NH3_emitted_kg;
        q.cap  = q.cap  + w * s.NH3_captured_kg;
        q.AS   = q.AS   + w * s.AS_kg;
        q.acid = q.acid + w * s.acid_kg;
        q.fan  = q.fan  + w * s.fan_kWh;
        q.air  = q.air  + w * s.air_treated_m3;
        L_min  = s.liquid_L_min;
        if ~isnan(s.dp_max_Pa), dp_max = max(dp_max, s.dp_max_Pa); end
    end
    N = p.flocks_per_year;
    eco.dp_max_Pa       = dp_max;                        % packed bed only (0 for spray)
    eco.location        = L.label;
    eco.scrubber_type   = p.scrubber_type;
    eco.NH3_emitted_kg  = N * q.emit;
    eco.NH3_captured_kg = N * q.cap;
    eco.capture_pct     = 100 * q.cap / q.emit;
    eco.AS_t            = N * q.AS / 1000;
    eco.acid_t          = N * q.acid / 1000;

    % --- utilities
    if isempty(p.run_days), run_days = p.flock_days; else, run_days = min(p.run_days, p.flock_days); end
    hours = N * run_days * 24;                           % scrubber runs only while switched on
    if is_packed
        P_pump = p.pump_head_packed;                     % MPa
        w_loss = p.water_loss_packed;
    else
        P_pump = p.P_nozzle;
        w_loss = p.water_loss_spray;
    end
    pump_kW      = (L_min/60000) * (P_pump*1e6) / p.pump_eff / 1000;
    eco.pump_kWh = pump_kW * hours;
    eco.fan_kWh  = N * q.fan;
    eco.water_m3 = N * q.air * w_loss * 1e-6;            % mL -> m3

    % --- costs, USD per year
    [C, Q_cap] = capital_cost(p);
    i = L.interest;  n = p.econ.life_years;
    CRF = i*(1+i)^n / ((1+i)^n - 1);
    eco.capacity_m3h = Q_cap;
    eco.capital      = C;
    c.capital     = C * CRF;
    c.maintenance = C * p.econ.maint_frac;
    c.acid        = eco.acid_t * L.price_acid;
    c.electricity = (eco.pump_kWh + eco.fan_kWh) * L.price_elec;
    c.water       = eco.water_m3 * L.price_water;
    c.labour      = p.econ.labour_h_week * (N*run_days/7) * L.wage;
    % liquid product: volume = AS / (AS concentration); AS g/L = N mol/L / 2 * M_AS
    AS_gL         = p.econ.product_N_molL / 2 * p.M_AS;
    eco.product_m3 = eco.AS_t * 1e6 / AS_gL / 1000;
    c.handling    = eco.product_m3 * p.econ.handling_USD_m3;
    eco.cost      = c;
    eco.total_cost = c.capital + c.maintenance + c.acid + c.electricity + c.water + c.labour + c.handling;
    eco.revenue    = eco.AS_t * L.price_AS * p.econ.AS_value_factor;
    eco.net_cost   = eco.total_cost - eco.revenue;

    % --- indicators
    birds = N * p.n_birds;
    eco.cost_per_kg_NH3    = eco.net_cost / eco.NH3_captured_kg;          % USD per kg NH3 removed
    eco.cost_per_bird      = eco.net_cost / birds;                        % USD per bird marketed
    eco.cost_per_kg_live   = eco.cost_per_bird / p.live_weight_kg;        % USD per kg live weight
    eco.breakeven_AS_price = eco.total_cost / eco.AS_t;                   % USD/t at which net cost = 0
    eco.N_recovered_g_bird = eco.NH3_captured_kg * 14.007/p.M_NH3 * 1000 / birds;  % g N per bird
    % Moore et al. 2025 states that the acid costs more than the N is worth. Per kg N captured:
    eco.acid_per_kgN   = p.M_H2SO4/(2*14.007) / 1000 * L.price_acid;                   % USD acid
    eco.AS_value_per_kgN = p.M_AS/(2*14.007) / 1000 * L.price_AS * p.econ.AS_value_factor;  % USD income
    eco.prices = L;
end

% =====================================================================
%  SENSITIVITY ANALYSIS (one input at a time)
% =====================================================================
function sen = run_sensitivity(p, loc)
% Each uncertain input is moved to its low and its high value while all the
% others stay at the base case. The output watched is the net cost per kg
% of NH3 removed. Inputs are then ranked by the size of the swing.
%
% This is a local, one-at-a-time analysis: it shows which inputs matter
% most, but not how they interact. (A global Sobol analysis comes later.)
    v = uncertain_inputs(p, loc);

    base = annual_economics(p, loc);
    n = size(v,1);
    y_lo = zeros(n,1);  y_hi = zeros(n,1);
    for k = 1:n
        q = setfield(p, v{k,2}{:}, v{k,3});   e = annual_economics(q, loc);  y_lo(k) = e.cost_per_kg_NH3;
        q = setfield(p, v{k,2}{:}, v{k,4});   e = annual_economics(q, loc);  y_hi(k) = e.cost_per_kg_NH3;
    end
    swing = abs(y_hi - y_lo);
    [~, order] = sort(swing, 'descend');

    sen.location = base.location;   sen.scrubber_type = p.scrubber_type;
    sen.base   = base.cost_per_kg_NH3;
    sen.label  = v(order,1);
    sen.low    = cell2mat(v(order,3));   sen.high = cell2mat(v(order,4));
    sen.y_low  = y_lo(order);            sen.y_high = y_hi(order);
    sen.swing  = swing(order);           sen.source = v(order,5);

    fprintf('\n=== SENSITIVITY, %s scrubber, %s: net cost per kg NH3 removed (base %.2f USD) ===\n', ...
            upper(p.scrubber_type), base.location, sen.base);
    fprintf('%-32s %9s %9s %11s %11s %8s   %s\n', 'Input', 'low', 'high', 'cost@low', 'cost@high', 'swing', 'range from');
    for k = 1:n
        fprintf('%-32s %9.3g %9.3g %11.2f %11.2f %8.2f   %s\n', sen.label{k}, sen.low(k), sen.high(k), ...
                sen.y_low(k), sen.y_high(k), sen.swing(k), sen.source{k});
    end
end

% =====================================================================
%  SIZING SWEEP
% =====================================================================
function siz = run_sizing(p, loc)
% Net cost and capture against scrubber size, everything else at base case.
% A small scrubber treats only the base ventilation (the concentrated air);
% a large one treats more of the summer air but costs more.
    if strcmpi(p.scrubber_type, 'packed')
        sizes = [1 1.5 2 3 4 6 9 12 18 24];          % m2 face area
    else
        sizes = [5 8 10 15 20 30 45 60 90 120];      % number of spray modules
    end
    n = numel(sizes);
    siz.location = '';  siz.scrubber_type = p.scrubber_type;  siz.size = sizes;
    siz.capacity_m3h = zeros(1,n);  siz.capture_pct = zeros(1,n);
    siz.cost_per_kg  = zeros(1,n);  siz.net_cost    = zeros(1,n);  siz.cost_per_bird = zeros(1,n);
    for k = 1:n
        q = p;
        if strcmpi(p.scrubber_type, 'packed'), q.bed_area = sizes(k); else, q.n_modules = sizes(k); end
        e = annual_economics(q, loc);
        siz.location        = e.location;
        siz.capacity_m3h(k) = e.capacity_m3h;
        siz.capture_pct(k)  = e.capture_pct;
        siz.cost_per_kg(k)  = e.cost_per_kg_NH3;
        siz.net_cost(k)     = e.net_cost;
        siz.cost_per_bird(k)= e.cost_per_bird;
    end
    [~, i] = min(siz.cost_per_kg);
    siz.best_index = i;

    fprintf('\n=== SIZING, %s scrubber, %s ===\n', upper(p.scrubber_type), siz.location);
    fprintf('%14s %12s %14s %16s %16s\n', 'Capacity m3/h', 'Capture %', 'Net cost USD/yr', 'USD per kg NH3', 'cents per bird');
    for k = 1:n
        if k == i, mark = '  <- lowest cost per kg'; else, mark = ''; end
        fprintf('%14.0f %12.1f %14.0f %16.2f %16.2f%s\n', siz.capacity_m3h(k), siz.capture_pct(k), ...
                siz.net_cost(k), siz.cost_per_kg(k), 100*siz.cost_per_bird(k), mark);
    end
end

function v = uncertain_inputs(p, loc)
% The uncertain inputs and their ranges, shared by the one-at-a-time and
% the global analysis.  {label, field path in p, low, high, source of range}
    L = location_data(loc);
    v = { ...
      'Fertiliser (AS) price, USD/t',   {'override','price_AS'},      L.range_AS(1),   L.range_AS(2),   '2026 market reports'; ...
      'Sulphuric acid price, USD/t',    {'override','price_acid'},    L.range_acid(1), L.range_acid(2), '2026 market reports'; ...
      'Electricity price, USD/kWh',     {'override','price_elec'},    L.range_elec(1), L.range_elec(2), 'tariffs'; ...
      'Interest rate',                  {'override','interest'},      0.04,            0.12,            'assumption'; ...
      'Emission slope, g/bird/d/d',     {'ER_slope'},                 0.024,           0.038,           'Wheeler et al. 2006, farm range'; ...
      'Ventilation scale',              {'Q_scale'},                  0.8,             1.2,             'assumption'; ...
      'Flocks per year',                {'flocks_per_year'},          5,               7,               'assumption'; ...
      'Liquid product value factor',    {'econ','AS_value_factor'},   0.5,             1.0,             'assumption'; ...
      'Capital cost factor',            {'econ','capex_factor'},      0.7,             p.capex_high,    'assumption; Moore et al. 2025'; ...
      'Pump efficiency',                {'pump_eff'},                 0.13,            0.60,            'Hadlocon et al. 2015a field ... good pump'; ...
      'Product handling, USD/m3',       {'econ','handling_USD_m3'},   0,               15,              'assumption'};
    if strcmpi(p.scrubber_type, 'packed')
        v(end+1,:) = {'Packing constant C_V', {'pack','C_V'}, 0.167, 0.35, 'Billet and Schultes 1999, Table 2'};
    else
        v(end+1,:) = {'Spray mass-transfer factor', {'kya_factor'}, 0.75, 1.36, 'Hadlocon et al. 2015c data; widened downward'};
    end
end

% =====================================================================
%  DESIGN OPTIMIZATION (exhaustive grid search)
% =====================================================================
function opt = run_optimization(p, loc)
% Every combination of the design variables is simulated (no optimization
% toolbox needed, and the global optimum on the grid is guaranteed).
%
%   Spray : number of modules, spray stages, nozzle pressure
%   Packed: face area, packing height, liquid load
%           constraint: packing pressure drop <= p.dp_limit
%
% Two questions are answered:
%   1. Which design has the lowest net cost per kg NH3 removed?
%   2. For a required capture (30 ... 90 %), which design has the lowest
%      net cost per year? (This is the cost-of-control curve.)
%
% Capital cost depends on capacity, number of spray stages and packing height; pH is not varied because at pH <= 5 the back-pressure
% is negligible and acid use is stoichiometric.
    is_packed = strcmpi(p.scrubber_type, 'packed');
    if is_packed
        g1 = [1 2 3 4 6 8 10 12 16 20 24];       names = {'area_m2', 'height_m', 'liquid_m3h_m2'};
        g2 = [0.3 0.45 0.6 0.9 1.2 1.5];
        g3 = [5 10 20];
    else
        g1 = [5 8 10 15 20 25 30 40 50 60 75 90 120];   names = {'modules', 'stages', 'P_MPa'};
        g2 = [1 2 3];
        g3 = [0.21 0.31 0.41 0.51 0.62];
    end
    n = numel(g1)*numel(g2)*numel(g3);
    D = zeros(n,3);  cap = zeros(n,1);  net = zeros(n,1);  cpk = zeros(n,1);
    cpb = zeros(n,1);  Qc = zeros(n,1);  dpm = zeros(n,1);
    k = 0;
    for a = g1
        for b = g2
            for c = g3
                k = k + 1;
                q = p;
                if is_packed
                    q.bed_area = a;  q.bed_height = b;  q.u_liquid = c;
                else
                    q.n_modules = a;  q.n_stages = b;  q.P_nozzle = c;
                end
                e = annual_economics(q, loc);
                D(k,:) = [a b c];
                cap(k) = e.capture_pct;       net(k) = e.net_cost;
                cpk(k) = e.cost_per_kg_NH3;   cpb(k) = e.cost_per_bird;
                Qc(k)  = e.capacity_m3h;      dpm(k) = e.dp_max_Pa;
                opt.location = e.location;
            end
        end
    end
    ok = dpm + p.dp_mist <= p.dp_limit;           % feasible: packing + mist eliminator within fan limit

    % 1. lowest cost per kg
    c1 = cpk;  c1(~ok) = inf;
    [~, ib] = min(c1);

    % 2. cost-of-control curve
    targets = 30:10:90;
    it = zeros(size(targets));
    for j = 1:numel(targets)
        c2 = net;  c2(~ok | cap < targets(j)) = inf;
        [m, i] = min(c2);
        if isinf(m), it(j) = 0; else, it(j) = i; end
    end

    opt.scrubber_type = p.scrubber_type;   opt.var_names = names;
    opt.design = D;  opt.capture_pct = cap;  opt.net_cost = net;  opt.cost_per_kg = cpk;
    opt.cost_per_bird = cpb;  opt.capacity_m3h = Qc;  opt.dp_max_Pa = dpm;  opt.feasible = ok;
    opt.best_index = ib;      opt.best_design = D(ib,:);
    opt.targets = targets;    opt.target_index = it;

    fprintf('\n=== OPTIMIZATION, %s scrubber, %s: %d designs simulated, %d feasible ===\n', ...
            upper(p.scrubber_type), opt.location, n, sum(ok));
    fprintf('%-22s %10s %10s %14s %12s %10s %12s %12s %10s\n', 'Goal', names{:}, 'Capacity m3/h', ...
            'Capture %', 'USD/kg NH3', 'Net USD/yr', 'cents/bird');
    row = @(lab, i) fprintf('%-22s %10.4g %10.4g %14.4g %12.0f %10.1f %12.2f %12.0f %10.2f\n', lab, ...
            D(i,1), D(i,2), D(i,3), Qc(i), cap(i), cpk(i), net(i), 100*cpb(i));
    row('Lowest cost per kg', ib);
    for j = 1:numel(targets)
        if it(j) > 0
            row(sprintf('Capture >= %d %%', targets(j)), it(j));
        else
            fprintf('%-22s %s\n', sprintf('Capture >= %d %%', targets(j)), '   not reachable within the grid and the pressure-drop limit');
        end
    end
end

% =====================================================================
%  GLOBAL UNCERTAINTY AND SENSITIVITY (Monte Carlo + Sobol indices)
% =====================================================================
function gsa = run_global(p, loc, N)
% All uncertain inputs are varied together, uniformly within their ranges.
%
%   Uncertainty : distribution of the net cost per kg NH3 removed, and the
%                 probability that the scrubber pays for itself (net cost < 0).
%   Sensitivity : Sobol indices by the Saltelli sampling scheme.
%                 S_i  (first order) = share of the output variance explained
%                                      by input i alone
%                 ST_i (total)       = share including its interactions
%                 Estimators (Saltelli et al. 2010; Jansen 1999):
%                   S_i  = mean( (fB - mean) .* (fABi - fA) ) / Var
%                   ST_i = mean( (fA - fABi).^2 ) / (2 Var)
%
% Sampling uses a scrambled Halton sequence (deterministic), so results are
% reproducible. Confidence intervals of the total indices
% by bootstrap over the N rows. Cost: N*(k+2) model runs, k = number of inputs.
    v = uncertain_inputs(p, loc);
    k = size(v,1);
    lo = cell2mat(v(:,3))';   hi = cell2mat(v(:,4))';

    H = halton_points(N, 2*k);
    A = repmat(lo, N, 1) + H(:, 1:k)     .* repmat(hi - lo, N, 1);
    B = repmat(lo, N, 1) + H(:, k+1:2*k) .* repmat(hi - lo, N, 1);

    fA = zeros(N,1);  fB = zeros(N,1);  fAB = zeros(N,k);
    for i = 1:N
        fA(i) = eval_sample(p, loc, v, A(i,:));
        fB(i) = eval_sample(p, loc, v, B(i,:));
    end
    for j = 1:k
        AB = A;  AB(:,j) = B(:,j);
        for i = 1:N
            fAB(i,j) = eval_sample(p, loc, v, AB(i,:));
        end
    end

    y   = [fA; fB];
    Var = var(y);
    y0  = mean(y);
    S   = zeros(k,1);  ST = zeros(k,1);
    for j = 1:k
        S(j)  = mean((fB - y0) .* (fAB(:,j) - fA)) / Var;   % centred: same expectation, less noise
        ST(j) = mean((fA - fAB(:,j)).^2) / (2*Var);
    end
    % bootstrap confidence intervals of the total indices (rows resampled, B = 500)
    B = 500;  STb = zeros(k, B);  st = 4242;  idx = zeros(N,1);
    for bi = 1:B
        for r = 1:N
            [st, u] = park_miller(st);  idx(r) = 1 + floor(u*N);
        end
        Vb = var([fA(idx); fB(idx)]);
        for j = 1:k
            STb(j,bi) = mean((fA(idx) - fAB(idx,j)).^2) / (2*Vb);
        end
    end
    STs = sort(STb, 2);
    ST_lo = STs(:, max(1, round(0.025*B)));   ST_hi = STs(:, round(0.975*B));
    [~, order] = sort(ST, 'descend');

    % distribution of the output: the 2N independent points of A and B only
    ys = sort(y);  m = numel(ys);
    pct = @(f) ys(max(1, min(m, round(f*m))));
    base = annual_economics(p, loc);

    gsa.location = base.location;   gsa.scrubber_type = p.scrubber_type;
    gsa.N = N;  gsa.runs = N*(k+2);
    gsa.y = y;                      % all Monte Carlo outputs, USD per kg NH3
    gsa.base = base.cost_per_kg_NH3;
    gsa.mean = mean(y);  gsa.p05 = pct(0.05);  gsa.p50 = pct(0.50);  gsa.p95 = pct(0.95);
    gsa.prob_profit = mean(y < 0);
    gsa.label = v(order,1);  gsa.S = S(order);  gsa.ST = ST(order);
    gsa.ST_lo = ST_lo(order);  gsa.ST_hi = ST_hi(order);

    fprintf('\n=== GLOBAL ANALYSIS, %s scrubber, %s: %d model runs ===\n', ...
            upper(p.scrubber_type), base.location, gsa.runs);
    fprintf('Net cost per kg NH3 removed (USD): base %.2f | mean %.2f | 5%% %.2f | median %.2f | 95%% %.2f\n', ...
            gsa.base, gsa.mean, gsa.p05, gsa.p50, gsa.p95);
    fprintf('Probability that the scrubber pays for itself (net cost < 0): %.1f %%\n', 100*gsa.prob_profit);
    fprintf('Distribution and probability from the %d independent sample points (matrices A and B).\n', 2*N);
    fprintf('%-32s %12s %12s %16s\n', 'Input', 'first order', 'total', 'total, 95% CI');
    for j = 1:k
        fprintf('%-32s %12.2f %12.2f %9.2f-%5.2f\n', gsa.label{j}, gsa.S(j), gsa.ST(j), gsa.ST_lo(j), gsa.ST_hi(j));
    end
    fprintf('Sum of first-order indices: %.2f (1.00 would mean no interactions).\n', sum(S));
end

function y = eval_sample(p, loc, v, x)
% Model output (USD per kg NH3 removed) for one set of uncertain inputs.
    q = p;
    for j = 1:size(v,1)
        q = setfield(q, v{j,2}{:}, x(j));
    end
    e = annual_economics(q, loc);
    y = e.cost_per_kg_NH3;
end

function H = halton_points(N, d)
% N points of the d-dimensional scrambled Halton sequence in (0,1), first 20
% skipped. Each dimension uses a fixed permutation of the digits 1..b-1
% (0 stays 0), drawn with the Park-Miller generator, so the scrambling is
% deterministic and the points are reproducible. Scrambling
% removes the correlation between dimensions of the plain Halton sequence.
    pr = [2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113];
    if d > numel(pr), error('Too many dimensions for halton_points.'); end
    H = zeros(N, d);
    for j = 1:d
        b = pr(j);
        perm = 0:b-1;                                    % perm(digit+1) = scrambled digit
        st = 1000 + 7919*j;
        for m = b:-1:3                                   % Fisher-Yates on digits 1..b-1
            [st, u] = park_miller(st);
            r = 2 + floor(u * (m - 1));                  % position 2..m
            tmp = perm(m);  perm(m) = perm(r);  perm(r) = tmp;
        end
        for i = 1:N
            n = i + 20;  f = 1;  r = 0;
            while n > 0
                f = f / b;
                r = r + f * perm(mod(n, b) + 1);
                n = floor(n / b);
            end
            H(i,j) = r;
        end
    end
end

function [st, u] = park_miller(st)
% Minimal-standard linear congruential generator; exact in double precision.
    st = mod(16807 * st, 2147483647);
    u  = st / 2147483647;
end

% =====================================================================
%  CHECK 4: BYPASS LOAD FRACTION AGAINST Melse et al. 2006 (broilers)
% =====================================================================
function bp = run_bypass_check(p)
% Melse et al. (2006) analysed hourly broiler-house data from the
% Netherlands: a scrubber that treats airflow only up to 50% of the maximum
% ventilation rate (1.75 m3/bird/h), with the rest bypassed, still receives
% about 85% of the NH3 load. Here the same quantity is computed from our
% daily house model, for each season:
%     load fraction = sum(E * min(Q, Q_set)/Q) / sum(E)
% Our model is daily, so short airflow peaks within a day are averaged out;
% this makes our fraction too high compared with an hourly analysis.
    day = (1:p.flock_days)';
    seasons = {'cold', 'mild', 'hot'};
    set_frac = [0.25 0.50 0.75];
    bp.season = seasons;  bp.set_frac = set_frac;  bp.load = zeros(3, numel(set_frac));
    bp.design = zeros(1,3);
    Q_design = p.n_modules * pi/4 * p.col_diameter^2 * p.v_max * 86400;   % m3/d, spray design
    fprintf('\n=== CHECK 4: share of NH3 load reaching a bypassed scrubber, against Melse et al. 2006 ===\n');
    fprintf('%-8s %12s %18s %18s %18s %18s\n', 'Season', 'Q max', 'set 25% of max', 'set 50% of max', 'set 75% of max', 'our design');
    fprintf('%-8s %12s %18s %18s %18s %18s\n', '', 'm3/bird/h', '', '(Melse et al. 2006: ~85%)', '', sprintf('%.2f m3/bird/h', Q_design/24/p.n_birds));
    for k = 1:3
        Q = ventilation_rate(day, seasons{k}, p);
        E = emission_rate(day, p) .* airflow_factor(day, seasons{k}, p);
        for j = 1:numel(set_frac)
            Qs = set_frac(j) * max(Q);
            bp.load(k,j) = 100 * sum(E .* min(Q, Qs) ./ Q) / sum(E);
        end
        bp.design(k) = 100 * sum(E .* min(Q, Q_design) ./ Q) / sum(E);
        fprintf('%-8s %12.2f %17.1f%% %17.1f%% %17.1f%% %17.1f%%\n', seasons{k}, max(Q)/24/p.n_birds, bp.load(k,:), bp.design(k));
    end
    fprintf(['Melse et al. 2006, broilers (Dutch data): a bypass set at 1.75 m3/bird/h (50%% of the hourly\n' ...
             'maximum, 3.5 m3/bird/h) still sends about 85%% of the NH3 load through the scrubber.\n' ...
             'Compare in absolute terms with the "our design" column (1.72 m3/bird/h): the Dutch\n' ...
             'climate lies between our cold and mild seasons. The "%% of max" columns are not\n' ...
             'directly comparable: our maximum is a daily mean, Melse et al. 2006 used hourly data, whose\n' ...
             'maximum is higher, so the same %% setting is a smaller scrubber in our model.\n']);
end

% =====================================================================
%  OPERATING STRATEGY: SCRUBBER ON ONLY IN THE FIRST PART OF THE FLOCK
% =====================================================================
function st = run_strategy(p, loc)
% Moore et al. 2025 suggests running a scrubber only during the first 3-4 weeks,
% when ventilation is low. Here the scrubber (same size) runs on days
% 1..run_days and is switched off afterwards (no pumping, no capture).
    days = [21 28 35 42 p.flock_days];
    n = numel(days);
    st.run_days = days;  st.capture_pct = zeros(1,n);  st.cost_per_kg = zeros(1,n);
    st.cost_per_bird = zeros(1,n);  st.net_cost = zeros(1,n);
    for k = 1:n
        q = p;  q.run_days = days(k);
        e = annual_economics(q, loc);
        st.location = e.location;
        st.capture_pct(k)   = e.capture_pct;
        st.cost_per_kg(k)   = e.cost_per_kg_NH3;
        st.cost_per_bird(k) = e.cost_per_bird;
        st.net_cost(k)      = e.net_cost;
    end
    fprintf('\n=== OPERATING STRATEGY, %s scrubber, %s: scrubber on from day 1 to day ... ===\n', upper(p.scrubber_type), st.location);
    fprintf('%10s %12s %16s %16s %16s\n', 'Last day', 'Capture %', 'Net cost USD/yr', 'USD per kg NH3', 'cents per bird');
    for k = 1:n
        fprintf('%10d %12.1f %16.0f %16.2f %16.2f\n', days(k), st.capture_pct(k), st.net_cost(k), st.cost_per_kg(k), 100*st.cost_per_bird(k));
    end
end

% =====================================================================
%  CAPITAL-COST SCENARIOS
% =====================================================================
function cx = run_capex(p, loc)
% Net cost when the capital cost is higher than our scaled estimate:
% factor 1 (base), 2, 3, and the European level of Moore et al. 2025 (250,000 USD
% per house for a full-capacity scrubber).
    C0 = capital_cost(p);
    f  = [1 2 3 250000/C0];
    lab = {'base estimate', '2 x base', '3 x base', 'European level (Moore et al. 2025)'};
    n = numel(f);
    cx.factor = f;  cx.label = lab;  cx.capital = zeros(1,n);  cx.cost_per_kg = zeros(1,n);
    cx.cost_per_bird = zeros(1,n);
    for k = 1:n
        q = p;  q.econ.capex_factor = f(k);
        e = annual_economics(q, loc);
        cx.location = e.location;
        cx.capital(k) = e.capital;  cx.cost_per_kg(k) = e.cost_per_kg_NH3;  cx.cost_per_bird(k) = e.cost_per_bird;
    end
    fprintf('\n=== CAPITAL-COST SCENARIOS, %s scrubber, %s ===\n', upper(p.scrubber_type), cx.location);
    fprintf('%-28s %8s %14s %16s %16s\n', 'Scenario', 'factor', 'Capital USD', 'USD per kg NH3', 'cents per bird');
    for k = 1:n
        fprintf('%-28s %8.1f %14.0f %16.2f %16.2f\n', lab{k}, f(k), cx.capital(k), cx.cost_per_kg(k), 100*cx.cost_per_bird(k));
    end
end

% =====================================================================
%  FURTHER SCENARIOS
% =====================================================================
function rv = run_further(p)
% 1. Climate against economics: net cost with the climate (season mix) of one
%    location and the prices of another.
% 2. Most and least favorable corner of the uncertainty ranges.
% 3. Does the optimal number of spray stages survive a cost per stage?
% 4. Day-night swing of the ventilation (hourly time step).
% 5. Emission slope calibrated to the measured emission of Moore et al. 2011.
    locs = {'us', 'gulf', 'developing'};  nm = {'US', 'Gulf', 'Developing'};
    q = p;  q.scrubber_type = 'spray';

    % ---- 1. climate x economics
    rv.cpk = zeros(3);  rv.cpb = zeros(3);
    for c = 1:3
        Lc = location_data(locs{c});
        for e = 1:3
            r = q;  r.override = struct('season_mix', Lc.season_mix);
            o = annual_economics(r, locs{e});
            rv.cpk(c,e) = o.cost_per_kg_NH3;  rv.cpb(c,e) = 100*o.cost_per_bird;
        end
    end
    fprintf('\n=== CLIMATE AGAINST ECONOMICS, SPRAY scrubber, base design ===\n');
    fprintf('Rows: climate (season mix) of ...; columns: prices of ...\n');
    fprintf('USD per kg NH3 removed (cents per bird in brackets)\n');
    fprintf('%-14s %18s %18s %18s\n', 'climate \\ prices', nm{:});
    for c = 1:3
        fprintf('%-14s', nm{c});
        for e = 1:3, fprintf('   %6.2f (%6.2f)', rv.cpk(c,e), rv.cpb(c,e)); end
        fprintf('\n');
    end
    fprintf('Climate effect with US prices: %.2f -> %.2f -> %.2f cents per bird.\n', rv.cpb(1,1), rv.cpb(2,1), rv.cpb(3,1));
    fprintf('Price effect with US climate:  %.2f -> %.2f -> %.2f cents per bird.\n', rv.cpb(1,1), rv.cpb(1,2), rv.cpb(1,3));

    % ---- 2. corners of the uncertainty ranges
    fprintf('\n=== MOST AND LEAST FAVORABLE COMBINATION OF ALL UNCERTAIN INPUTS, SPRAY ===\n');
    fprintf('%-14s %14s %14s %14s %16s\n', 'Location', 'base USD/kg', 'best USD/kg', 'worst USD/kg', 'best net USD/yr');
    rv.corner = zeros(3,3);
    for e = 1:3
        v = uncertain_inputs(q, locs{e});  k = size(v,1);
        lo = cell2mat(v(:,3))';  hi = cell2mat(v(:,4))';
        b0 = annual_economics(q, locs{e});
        xb = zeros(1,k);  xw = zeros(1,k);
        x0 = zeros(1,k);
        L = location_data(locs{e});
        for j = 1:k
            if strcmp(v{j,2}{1}, 'override')
                x0(j) = L.(v{j,2}{2});     % prices come from the location data
            else
                x0(j) = getfield(q, v{j,2}{:}); %#ok<GFLD>
            end
        end
        for j = 1:k
            xl = x0;  xl(j) = lo(j);  xh = x0;  xh(j) = hi(j);
            yl = eval_sample(q, locs{e}, v, xl);  yh = eval_sample(q, locs{e}, v, xh);
            if yl <= yh, xb(j) = lo(j); xw(j) = hi(j); else, xb(j) = hi(j); xw(j) = lo(j); end
        end
        yb = eval_sample(q, locs{e}, v, xb);  yw = eval_sample(q, locs{e}, v, xw);
        r = q;  for j = 1:k, r = setfield(r, v{j,2}{:}, xb(j)); end
        ob = annual_economics(r, locs{e});
        rv.corner(e,:) = [b0.cost_per_kg_NH3 yb yw];
        fprintf('%-14s %14.2f %14.2f %14.2f %16.0f\n', nm{e}, b0.cost_per_kg_NH3, yb, yw, ob.net_cost);
    end
    fprintf('"best" = every input at the end of its range that lowers the cost (one at a time).\n');

    % ---- 3. cost per spray stage and the optimal design
    fprintf('\n=== OPTIMAL SPRAY DESIGN AGAINST THE CAPITAL COST OF A SPRAY STAGE ===\n');
    fprintf('%-12s %10s | %-28s | %-28s\n', 'Location', 'stage cost', 'lowest cost per kg', 'cheapest design for >= 70 %');
    g1 = [5 8 10 15 20 25 30 40 50 60 75 90 120];  g3 = [0.21 0.31 0.41 0.51 0.62];
    sc = [0 0.15 0.25];
    rv.stage = cell(3, numel(sc));
    for e = 1:3
        for t = 1:numel(sc)
            best = [inf 0 0 0 0];  b70 = [inf 0 0 0 0];
            for a = g1
                for b = 1:3
                    for c = g3
                        r = q;  r.n_modules = a;  r.n_stages = b;  r.P_nozzle = c;  r.econ.stage_cost = sc(t);
                        o = annual_economics(r, locs{e});
                        if o.cost_per_kg_NH3 < best(1), best = [o.cost_per_kg_NH3 a b c o.capture_pct]; end
                        if o.capture_pct >= 70 && o.net_cost < b70(1), b70 = [o.net_cost a b c 100*o.cost_per_bird]; end
                    end
                end
            end
            rv.stage{e,t} = [best; b70];
            fprintf('%-12s %10.2f | %3d mod, %d st, %.2f MPa, %5.2f | %3d mod, %d st, %.2f MPa, %5.2f c\n', ...
                nm{e}, sc(t), best(2), best(3), best(4), best(1), b70(2), b70(3), b70(4), b70(5));
        end
    end
    fprintf('Stage cost = fraction of capital added (removed) per spray stage above (below) three.\n');

    % ---- 4. day-night swing of the ventilation
    fprintf('\n=== DAY-NIGHT SWING OF THE VENTILATION (hourly time step), SPRAY, base design ===\n');
    fprintf('%-12s %10s %12s %14s %14s\n', 'Location', 'amplitude', 'capture %', 'USD per kg', 'cents per bird');
    amps = [0 0.3 0.5];
    rv.diurnal = zeros(3, numel(amps), 3);
    for e = 1:3
        for t = 1:numel(amps)
            r = q;  r.diurnal_amp = amps(t);
            o = annual_economics(r, locs{e});
            rv.diurnal(e,t,:) = [o.capture_pct o.cost_per_kg_NH3 100*o.cost_per_bird];
            fprintf('%-12s %10.1f %12.1f %14.2f %14.2f\n', nm{e}, amps(t), o.capture_pct, o.cost_per_kg_NH3, 100*o.cost_per_bird);
        end
    end
    fprintf('Amplitude A: hourly airflow = daily mean x (1 + A cos(2 pi (h - 15)/24)).\n');

    % ---- 5. emission slope calibrated to Moore et al. 2011
    E_meas = mean([13.6 16.3 14.0]);
    day = (1:p.flock_days)';  E_mod = mean(emission_rate(day, p));
    slope = p.ER_slope * E_meas / E_mod;
    fprintf('\n=== EMISSION SLOPE CALIBRATED TO Moore et al. 2011 (%.4f instead of %.4f g/bird/d/d) ===\n', slope, p.ER_slope);
    fprintf('%-12s %14s %14s %14s %14s\n', 'Location', 'USD/kg base', 'USD/kg calib.', 'c/bird base', 'c/bird calib.');
    rv.calib = zeros(3,4);
    for e = 1:3
        o0 = annual_economics(q, locs{e});
        r = q;  r.ER_slope = slope;  o1 = annual_economics(r, locs{e});
        rv.calib(e,:) = [o0.cost_per_kg_NH3 o1.cost_per_kg_NH3 100*o0.cost_per_bird 100*o1.cost_per_bird];
        fprintf('%-12s %14.2f %14.2f %14.2f %14.2f\n', nm{e}, rv.calib(e,:));
    end
end

% =====================================================================
%  JOURNAL FIGURES
% =====================================================================
% Each figure is drawn at its printed size in cm (Elsevier: 9 cm single column,
% 14 cm 1.5 columns, 19 cm full width) with 8 pt Arial text. Panels are placed in
% cm, labelled (a), (b), ..., and the same colors are used for the same thing in
% every figure. Saved as vector PDF (journal), EMF (Word, Windows), PNG 600 dpi,
% and each panel as its own file in the subfolder 'panels'.
function make_figures(D, save_flag, outdir)
    if nargin < 2, save_flag = false; end
    if nargin < 3, outdir = fullfile(pwd, '..', 'figures'); end
    S = fig_style();
    F = {};                                   % {handle, file name, {axes, panel names}}
    F(end+1,:) = fig_checks(D, S);
    F(end+1,:) = fig_flock(D, S, 1, 'Fig3_flock_spray');
    F(end+1,:) = fig_cost(D, S);
    F(end+1,:) = fig_control(D, S);
    F(end+1,:) = fig_tornado(D, S);
    if isfield(D, 'gsa')
        F(end+1,:) = fig_distribution(D, S);
        F(end+1,:) = fig_sobol(D, S);
    end
    F(end+1,:) = fig_scenarios(D, S, isfield(D, 'gsa'));
    F(end+1,:) = fig_sizing(D, S);
    F(end+1,:) = fig_flock(D, S, 2, 'FigS2_flock_packed');
    if save_flag
        if ~exist(outdir, 'dir'), mkdir(outdir); end
        pdir = fullfile(outdir, 'panels');
        if ~exist(pdir, 'dir'), mkdir(pdir); end
        for k = 1:size(F, 1)
            save_figure(F{k,1}, fullfile(outdir, F{k,2}));
            P = F{k,3};
            for j = 1:size(P, 1)
                save_panel(P{j,1}, fullfile(pdir, [F{k,2} '_' P{j,2}]), S);
            end
        end
        fprintf('\nFigures saved in %s (PDF vector, EMF, PNG 600 dpi; single panels in panels\\)\n', outdir);
    end
end

function g = hzeros(n)
% Array for graphics handles: graphics objects (gobjects).
    try, g = gobjects(1, n); catch, g = zeros(1, n); end
end

function S = fig_style()
    S.font = 'Arial';  S.fs = 8;  S.lw = 1.2;  S.axlw = 0.6;  S.ms = 4;
    S.loc  = [0.00 0.30 0.60;      % US          dark blue
              0.85 0.40 0.05;      % Gulf        orange
              0.20 0.60 0.30];     % Developing  green
    S.locname = {'US', 'Gulf', 'Developing'};
    S.locmark = {'o', 's', '^'};
    S.sea  = [0.10 0.45 0.80;      % cold  blue
              0.40 0.65 0.15;      % mild  green
              0.85 0.20 0.15];     % hot   red
    S.cost = [0.20 0.40 0.65;      % capital
              0.60 0.75 0.90;      % maintenance
              0.95 0.65 0.15;      % acid
              0.55 0.35 0.65;      % electricity
              0.35 0.75 0.80;      % water
              0.70 0.80 0.40;      % labor
              0.75 0.30 0.30];     % handling
    S.bar  = [0.35 0.55 0.80];
end

function h = new_fig(name, W, H)
% Figure window of W x H cm; the saved page has exactly this size.
    h = figure('Name', name, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 W H], ...
               'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual', ...
               'PaperPosition', [0 0 W H], 'PaperSize', [W H], 'InvertHardcopy', 'off');
    try, set(h, 'Renderer', 'painters'); catch, end
end

function ax = cm_axes(h, x, y, w, ht, S)
% Axes placed in cm from the lower-left corner of the figure.
    pos = get(h, 'Position');  W = pos(3);  H = pos(4);
    ax = axes('Parent', h, 'Units', 'normalized', 'Position', [x/W y/H w/W ht/H]);
    hold(ax, 'on');  box(ax, 'off');   % axes on the left and bottom only (journal style)
    set(ax, 'FontName', S.font, 'FontSize', S.fs, 'LineWidth', S.axlw, 'TickDir', 'out', ...
            'TickLength', [0.015 0.015], 'Layer', 'top');
end

function plabel(ax, letter, txt, S)
% Panel label "(a) text" above the top-left corner of the axes.
    if nargin < 3, txt = ''; end
    t = ['(' letter ')'];
    if ~isempty(txt), t = [t ' ' txt]; end
    text('Parent', ax, 'Units', 'normalized', 'Position', [0 1.03], 'String', t, ...
         'FontName', S.font, 'FontSize', S.fs + 1, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

function style_labels(ax, xl, yl, S)
    if ~isempty(xl), xlabel(ax, xl, 'FontName', S.font, 'FontSize', S.fs); end
    if ~isempty(yl), ylabel(ax, yl, 'FontName', S.font, 'FontSize', S.fs); end
end

function lg = top_legend(h, hs, names, S, ycm, ncol_hint)
% One shared horizontal legend centered at height ycm (cm) of the figure.
    lg = legend(hs, names, 'Orientation', 'horizontal', 'FontName', S.font, 'FontSize', S.fs);
    set(lg, 'Box', 'off', 'Units', 'centimeters');
    pos = get(h, 'Position');  lp = get(lg, 'Position');
    set(lg, 'Position', [(pos(3) - lp(3))/2, ycm, lp(3), lp(4)]);
end

% ---------------------------------------------------------------- Fig. 2
function out = fig_checks(D, S)
    h = new_fig('Fig. 2 Model checks', 14, 7.4);
    val = D.val;  pk = D.pk;
    a1 = cm_axes(h, 1.4, 1.2, 5.3, 5.3, S);
    isEta = val.measured <= 100 & val.model <= 100;
    in = isEta & val.in_range;  ex = isEta & ~val.in_range;
    g1 = plot(a1, [30 100], [30 100], 'k-', 'LineWidth', 0.8);
    g2 = plot(a1, [30 100], [40 110], 'k:', 'LineWidth', 0.8);  plot(a1, [40 100], [30 90], 'k:', 'LineWidth', 0.8);
    m1 = plot(a1, val.measured(in), val.model(in), 'o', 'Color', S.loc(1,:), 'MarkerFaceColor', S.loc(1,:), 'MarkerSize', 5);
    m2 = plot(a1, val.measured(ex), val.model(ex), 's', 'Color', S.loc(2,:), 'LineWidth', 1, 'MarkerSize', 5);
    axis(a1, [40 100 30 100]);  style_labels(a1, 'Measured removal (%)', 'Model removal (%)', S);   % y from 30: room for the legend
    lg = legend(a1, [m1 m2 g1 g2], {'in range', 'extrapolated', '1:1', '\pm10 points'}, ...
           'Location', 'southeast', 'FontSize', S.fs - 1, 'Box', 'off');
    try, set(lg, 'ItemTokenSize', [14 9]); catch, end     % shorter line samples (MATLAB)
    plabel(a1, 'a', 'Spray scrubber', S);

    a2 = cm_axes(h, 8.2, 1.2, 5.3, 5.3, S);
    plot(a2, [40 100], [40 100], 'k-', 'LineWidth', 0.8);
    plot(a2, [40 100], [50 110], 'k:', 'LineWidth', 0.8);  plot(a2, [40 100], [30 90], 'k:', 'LineWidth', 0.8);
    n1 = plot(a2, pk.measured, pk.model(:,1), 'o', 'Color', S.loc(1,:), 'MarkerFaceColor', S.loc(1,:), 'MarkerSize', 5);
    n2 = plot(a2, pk.measured, pk.model(:,2), 's', 'Color', S.loc(2,:), 'LineWidth', 1, 'MarkerSize', 5);
    axis(a2, [40 100 30 100]);  style_labels(a2, 'Measured removal (%)', 'Model removal (%)', S);
    legend(a2, [n1 n2], {'C_V = 0.167', 'C_V = 0.35'}, 'Location', 'southeast', 'FontSize', S.fs - 1, 'Box', 'off');
    plabel(a2, 'b', 'Packed bed', S);
    out = {h, 'Fig2_model_checks', {a1, 'a_spray'; a2, 'b_packed'}};
end

% ---------------------------------------------------------------- Fig. 3 and S2
function out = fig_flock(D, S, t, fname)
    nm = {'Fig. 3 One flock, spray scrubber', 'Fig. S2 One flock, packed bed'};
    h = new_fig(nm{t}, 19, 13.4);
    res = D.res;  ns = numel(D.seasons);
    X = [1.5 11.0];  Y = [7.3 1.1];  w = 7.2;  ht = 4.6;
    ax = hzeros(4);
    for i = 1:4
        ax(i) = cm_axes(h, X(mod(i-1,2)+1), Y(ceil(i/2)), w, ht, S);
    end
    hs = hzeros(ns);
    for k = 1:ns
        r = res{t,k};  c = S.sea(k,:);
        hs(k) = plot(ax(1), r.day, r.Q_m3h/1000, '-', 'Color', c, 'LineWidth', S.lw);
        plot(ax(2), r.day, r.C_in, '-', 'Color', c, 'LineWidth', S.lw);
        plot(ax(3), r.day, 100*r.capture, '-', 'Color', c, 'LineWidth', S.lw);
        plot(ax(4), r.day, r.AS, '-', 'Color', c, 'LineWidth', S.lw);
    end
    r = res{t,1};
    hc = plot(ax(1), r.day, r.summary.Q_cap_m3h/1000*ones(size(r.day)), 'k--', 'LineWidth', 1);
    plot(ax(4), r.day, r.E*132.14/(2*17.031), 'k:', 'LineWidth', 1.2);
    hp = plot(ax(1), NaN, NaN, 'k:', 'LineWidth', 1.2);   % legend entry for panel (d); all legend items must sit in one axes
    ylim(ax(3), [0 100]);
    for i = 1:4, xlim(ax(i), [0 50]); end
    style_labels(ax(1), '', 'Exhaust airflow (1000 m^3 h^{-1})', S);
    style_labels(ax(2), '', 'NH_3 at scrubber inlet (ppmv)', S);
    style_labels(ax(3), 'Bird age (d)', 'NH_3 captured (%)', S);
    style_labels(ax(4), 'Bird age (d)', 'Ammonium sulphate (kg d^{-1})', S);
    plabel(ax(1), 'a', 'House airflow', S);   plabel(ax(2), 'b', 'Inlet concentration', S);
    plabel(ax(3), 'c', 'Capture', S);         plabel(ax(4), 'd', 'Product', S);
    top_legend(h, [hs hc hp], [{'cold season', 'mild season', 'hot season'}, ...
               {'scrubber capacity', 'all NH_3 captured'}], S, 12.7);
    out = {h, fname, {ax(1), 'a_airflow'; ax(2), 'b_inlet'; ax(3), 'c_capture'; ax(4), 'd_product'}};
end

% ---------------------------------------------------------------- Fig. 4
function out = fig_cost(D, S)
    h = new_fig('Fig. 4 Annual cost and income', 19, 7.8);
    nl = numel(D.locations);  tt = {'Spray scrubber', 'Packed bed'};  let = 'ab';
    ax = hzeros(2);
    for t = 1:2
        ax(t) = cm_axes(h, 1.4 + (t-1)*6.9, 1.2, 5.6, 5.8, S);
        M = zeros(nl, 7);  R = zeros(nl, 1);
        for k = 1:nl
            c = D.eco{t,k}.cost;
            M(k,:) = [c.capital c.maintenance c.acid c.electricity c.water c.labour c.handling];
            R(k) = D.eco{t,k}.revenue;
        end
        hb = bar(ax(t), 1:nl, M/1000, 0.6, 'stacked');
        for j = 1:numel(hb), set(hb(j), 'FaceColor', S.cost(j,:), 'EdgeColor', 'k', 'LineWidth', 0.4); end
        hd = plot(ax(t), 1:nl, R/1000, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
        set(ax(t), 'XTick', 1:nl, 'XTickLabel', S.locname);  xlim(ax(t), [0.4 nl+0.6]);
        style_labels(ax(t), '', 'Annual cost or income (1000 USD)', S);
        plabel(ax(t), let(t), tt{t}, S);
    end
    yl = max([get(ax(1), 'YLim') get(ax(2), 'YLim')]);
    for t = 1:2, ylim(ax(t), [0 yl]); end            % same scale for both panels
    lg = legend(ax(2), [reshape(hb(end:-1:1), 1, []) hd], {'handling', 'labour', 'water', 'electricity', 'acid', ...
                'maintenance', 'capital', 'fertiliser income'}, 'FontName', S.font, 'FontSize', S.fs);
    set(lg, 'Box', 'off', 'Units', 'centimeters');  lp = get(lg, 'Position');
    set(lg, 'Position', [15.0, 1.2 + 5.8 - lp(4), lp(3), lp(4)]);
    out = {h, 'Fig4_annual_cost', {ax(1), 'a_spray'; ax(2), 'b_packed'}};
end

% ---------------------------------------------------------------- Fig. 5
function out = fig_control(D, S)
    h = new_fig('Fig. 5 Cost of ammonia control', 19, 7.6);
    tt = {'Spray scrubber', 'Packed bed'};  let = 'ab';  ax = hzeros(2);  hs = hzeros(3);
    for t = 1:2
        ax(t) = cm_axes(h, 1.4 + (t-1)*9.5, 1.1, 7.6, 5.2, S);
        for k = 1:numel(D.locations)
            o = D.opt{t,k};  ix = o.target_index(o.target_index > 0);
            hs(k) = plot(ax(t), o.capture_pct(ix), 100*o.cost_per_bird(ix), ['-' S.locmark{k}], ...
                         'Color', S.loc(k,:), 'MarkerFaceColor', S.loc(k,:), 'LineWidth', S.lw, 'MarkerSize', S.ms);
        end
        xlim(ax(t), [25 95]);
        style_labels(ax(t), 'NH_3 captured (% of house emission)', 'Net cost (US cents per bird)', S);
        plabel(ax(t), let(t), tt{t}, S);
    end
    top_legend(h, hs, S.locname, S, 6.9);
    out = {h, 'Fig5_cost_of_control', {ax(1), 'a_spray'; ax(2), 'b_packed'}};
end

% ---------------------------------------------------------------- Fig. 6
function out = fig_tornado(D, S)
    sen = D.sen;  nl = numel(sen);
    lab = sen{1}.label;  n = numel(lab);  ms = zeros(n,1);
    for k = 1:nl
        for j = 1:n, ms(j) = ms(j) + sen{k}.swing(strcmp(sen{k}.label, lab{j})); end
    end
    [~, o] = sort(ms, 'descend');  lab = lab(o);
    h = new_fig('Fig. 6 One-at-a-time sensitivity', 19, 9.6);
    ax = hzeros(nl);  let = 'abc';
    for k = 1:nl
        ax(k) = cm_axes(h, 4.4 + (k-1)*4.95, 1.1, 4.5, 7.3, S);
        s = sen{k};
        for j = 1:n
            i = find(strcmp(s.label, lab{j}));  y = n - j + 1;
            lo = min(s.y_low(i), s.y_high(i));  hi = max(s.y_low(i), s.y_high(i));
            patch([lo hi hi lo], [y-0.3 y-0.3 y+0.3 y+0.3], S.bar, 'Parent', ax(k), 'EdgeColor', 'none');
            h1 = plot(ax(k), s.y_low(i),  y, 'v', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4.5);
            h2 = plot(ax(k), s.y_high(i), y, '^', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 4.5);
        end
        h3 = plot(ax(k), [s.base s.base], [0.3 n+0.7], 'k--', 'LineWidth', 0.8);
        ylim(ax(k), [0.3 n+0.7]);
        xr = [min([s.y_low(:); s.y_high(:)]) max([s.y_low(:); s.y_high(:)])];
        xlim(ax(k), [floor(xr(1) - 0.05*diff(xr)) ceil(xr(2) + 0.05*diff(xr))]);
        if k == 1
            lab_print = strrep(strrep(strrep(lab, 'g/bird/d/d', 'g bird^{-1} d^{-2}'), 'USD/t', 'USD t^{-1}'), 'USD/kWh', 'USD kWh^{-1}');
            lab_print = strrep(lab_print, 'USD/m3', 'USD m^{-3}');
            set(ax(k), 'YTick', 1:n, 'YTickLabel', lab_print(end:-1:1));
        else
            set(ax(k), 'YTick', 1:n, 'YTickLabel', {});
        end
        style_labels(ax(k), 'USD per kg NH_3 removed', '', S);
        plabel(ax(k), let(k), S.locname{k}, S);
    end
    top_legend(h, [h1 h2 h3], {'input at low value', 'input at high value', 'base case'}, S, 9.0);
    out = {h, 'Fig6_sensitivity_tornado', {ax(1), 'a_US'; ax(2), 'b_Gulf'; ax(3), 'c_Developing'}};
end

% ---------------------------------------------------------------- Fig. 7
function out = fig_distribution(D, S)
    g = D.gsa;  nl = numel(g);
    h = new_fig('Fig. 7 Distribution of the net cost', 19, 6.9);
    ax = hzeros(nl);  let = 'abc';  ymax = 0;
    for k = 1:nl
        ax(k) = cm_axes(h, 1.3 + (k-1)*6.05, 1.1, 5.2, 4.3, S);
        y = g{k}.y;  edges = floor(min(y)):0.5:ceil(max(y));
        c = histc(y, edges);  c = c(:)';  c = c(1:end-1);  ctr = edges(1:end-1) + 0.25;
        bar(ax(k), ctr, 100*c/numel(y), 1, 'FaceColor', S.loc(k,:), 'EdgeColor', 'w', 'LineWidth', 0.3);
        ymax = max(ymax, max(100*c/numel(y)));
        xlim(ax(k), [edges(1) edges(end)]);
        style_labels(ax(k), 'USD per kg NH_3 removed', '', S);
        if k == 1, style_labels(ax(k), '', 'Share of samples (%)', S); end
        plabel(ax(k), let(k), S.locname{k}, S);
    end
    for k = 1:nl
        ylim(ax(k), [0 1.25*ymax]);
        hb = plot(ax(k), [g{k}.base g{k}.base], [0 1.25*ymax], 'k--', 'LineWidth', 1);
        hm = plot(ax(k), [g{k}.p50 g{k}.p50],   [0 1.25*ymax], 'k-',  'LineWidth', 1);
        xl = get(ax(k), 'XLim');
        if k > 1, set(ax(k), 'YTickLabel', {}); end     % same y scale: numbers on panel (a) only
        text('Parent', ax(k), 'Position', [xl(2) - 0.03*diff(xl), 1.19*ymax], ...
             'String', {sprintf('median %.2f', g{k}.p50), sprintf('5-95%%: %.2f-%.2f', g{k}.p05, g{k}.p95)}, ...
             'FontName', S.font, 'FontSize', S.fs - 1, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end
    top_legend(h, [hb hm], {'base case', 'median'}, S, 6.2);
    out = {h, 'Fig7_cost_distribution', {ax(1), 'a_US'; ax(2), 'b_Gulf'; ax(3), 'c_Developing'}};
end

% ---------------------------------------------------------------- Fig. 8
function out = fig_sobol(D, S)
    g = D.gsa;  nl = numel(g);
    lab = g{1}.label;  n = numel(lab);  ST = zeros(n, nl);  LO = ST;  HI = ST;
    for k = 1:nl
        for j = 1:n
            i = strcmp(g{k}.label, lab{j});
            ST(j,k) = max(g{k}.ST(i), 0);
            if isfield(g{k}, 'ST_lo'), LO(j,k) = max(g{k}.ST_lo(i), 0); HI(j,k) = g{k}.ST_hi(i);
            else, LO(j,k) = ST(j,k); HI(j,k) = ST(j,k); end
        end
    end
    [~, o] = sort(mean(ST, 2), 'descend');  lab = lab(o);  ST = ST(o,:);  LO = LO(o,:);  HI = HI(o,:);
    h = new_fig('Fig. 8 Total Sobol indices', 14, 9.9);
    ax = cm_axes(h, 4.6, 1.1, 9.0, 8.4, S);
    bh = 0.26;  hs = hzeros(nl);
    for j = 1:n
        yc = n - j + 1;
        for k = 1:nl
            y = yc + (2 - k)*bh;                      % US on top within each group
            hs(k) = patch([0 ST(j,k) ST(j,k) 0], [y-bh/2 y-bh/2 y+bh/2 y+bh/2], S.loc(k,:), 'Parent', ax, 'EdgeColor', 'none');
            plot(ax, [LO(j,k) HI(j,k)], [y y], 'k-', 'LineWidth', 0.6);
        end
    end
    ylim(ax, [0.4 n+0.6]);  xlim(ax, [0 ceil(10*max(HI(:)))/10]);
    lab_print = strrep(strrep(strrep(lab, 'g/bird/d/d', 'g bird^{-1} d^{-2}'), 'USD/t', 'USD t^{-1}'), 'USD/kWh', 'USD kWh^{-1}');
    lab_print = strrep(lab_print, 'USD/m3', 'USD m^{-3}');
    set(ax, 'YTick', 1:n, 'YTickLabel', lab_print(end:-1:1));
    style_labels(ax, 'Total Sobol index (share of the variance of the net cost)', '', S);
    lg = legend(ax, hs, S.locname, 'Location', 'southeast', 'FontName', S.font, 'FontSize', S.fs);
    set(lg, 'Box', 'off');
    out = {h, 'Fig8_sobol_indices', {}};
end

% ---------------------------------------------------------------- Fig. 9
function out = fig_scenarios(D, S, has_gsa)
    num = 9;  if ~has_gsa, num = 7; end
    h = new_fig(sprintf('Fig. %d Scenarios', num), 19, 7.6);
    ax = hzeros(2);  hs = hzeros(3);
    ax(1) = cm_axes(h, 1.4, 1.1, 7.6, 5.2, S);
    ax(2) = cm_axes(h, 10.9, 1.1, 7.6, 5.2, S);
    for k = 1:numel(D.locations)
        c = S.loc(k,:);  mk = ['-' S.locmark{k}];
        hs(k) = plot(ax(1), D.cx{k}.capital/1000, D.cx{k}.cost_per_kg, mk, 'Color', c, 'MarkerFaceColor', c, 'LineWidth', S.lw, 'MarkerSize', S.ms);
        plot(ax(2), D.st{k}.run_days, D.st{k}.cost_per_kg, mk, 'Color', c, 'MarkerFaceColor', c, 'LineWidth', S.lw, 'MarkerSize', S.ms);
    end
    style_labels(ax(1), 'Capital cost (1000 USD)', 'Net cost (USD per kg NH_3 removed)', S);
    style_labels(ax(2), 'Scrubber switched off after day', 'Net cost (USD per kg NH_3 removed)', S);
    plabel(ax(1), 'a', 'Capital cost up to the European level', S);
    plabel(ax(2), 'b', 'Running for part of the flock', S);
    top_legend(h, hs, S.locname, S, 6.9);
    out = {h, sprintf('Fig%d_scenarios', num), {ax(1), 'a_capital'; ax(2), 'b_part_flock'}};
end

% ---------------------------------------------------------------- Fig. S1
function out = fig_sizing(D, S)
    h = new_fig('Fig. S1 Sizing', 19, 13.4);
    X = [1.5 11.0];  Y = [7.3 1.1];  w = 7.2;  ht = 4.6;  ax = hzeros(4);  hs = hzeros(3);
    tt = {'Spray scrubber: cost', 'Packed bed: cost', 'Spray scrubber: capture', 'Packed bed: capture'};
    for t = 1:2
        ax(t)   = cm_axes(h, X(t), Y(1), w, ht, S);
        ax(t+2) = cm_axes(h, X(t), Y(2), w, ht, S);
        for k = 1:numel(D.locations)
            z = D.siz{t,k};  c = S.loc(k,:);  mk = ['-' S.locmark{k}];
            hs(k) = plot(ax(t), z.capacity_m3h/1000, z.cost_per_kg, mk, 'Color', c, 'MarkerFaceColor', c, 'LineWidth', S.lw, 'MarkerSize', S.ms);
            plot(ax(t+2), z.capacity_m3h/1000, z.capture_pct, mk, 'Color', c, 'MarkerFaceColor', c, 'LineWidth', S.lw, 'MarkerSize', S.ms);
        end
        style_labels(ax(t), '', 'USD per kg NH_3 removed', S);
        style_labels(ax(t+2), 'Scrubber capacity (1000 m^3 h^{-1})', 'NH_3 captured (%)', S);
        ylim(ax(t+2), [0 100]);
    end
    let = 'abcd';
    for i = 1:4, plabel(ax(i), let(i), tt{i}, S); end
    top_legend(h, hs, S.locname, S, 12.7);
    out = {h, 'FigS1_sizing', {ax(1), 'a'; ax(2), 'b'; ax(3), 'c'; ax(4), 'd'}};
end

% ---------------------------------------------------------------- saving
function save_figure(h, base)
    is_oct = exist('OCTAVE_VERSION', 'builtin') > 0;
    try                                            % MATLAB: no axes toolbar in the file
        a = findall(h, 'Type', 'axes');
        for j = 1:numel(a), a(j).Toolbar.Visible = 'off'; end
    catch
    end
    drawnow;
    try
        if is_oct
            print(h, [base '.pdf'], '-dpdf');
        else
            print(h, [base '.pdf'], '-dpdf', '-vector');               % journal (vector)
            if ispc, print(h, [base '.emf'], '-dmeta', '-vector'); end % Word (vector)
        end
        print(h, [base '.png'], '-dpng', '-r600');
    catch err
        fprintf('Could not save %s: %s\n', base, err.message);
    end
end

function save_panel(ax, base, S)
% One panel on its own page: the axes (with its own legend, if any) copied to a
% new figure with margins for the labels.
    pos = get(ax, 'Position');  hp = get(ax, 'Parent');  set(hp, 'Units', 'centimeters');
    fp = get(hp, 'Position');  w = pos(3)*fp(3);  ht = pos(4)*fp(4);
    W = w + 2.4;  H = ht + 2.0;
    h = new_fig('panel', W, H);
    lg = [];
    try, lg = get(ax, 'Legend'); catch, end
    try
        if ~isempty(lg) && isgraphics(lg) && ~strcmp(get(lg, 'Location'), 'none')
            c = copyobj([lg ax], h);  a2 = c(2);
        else
            a2 = copyobj(ax, h);
        end
        set(a2, 'Units', 'centimeters', 'Position', [1.8 1.2 w ht]);
        save_figure(h, base);
    catch err
        fprintf('Could not save panel %s: %s\n', base, err.message);
    end
    close(h);
end
