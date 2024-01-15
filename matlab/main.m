% Use this single script to generate synthetic DCE-MRI data:
%   - save input data (table 'ROI_param.csv', indexed phantom 'phantom.png' 
%       and coil sensitivities 'sensivitities.mat') to \data folder
%   - set simulation and acquisition parameters at the beggining of this 
%       script and run to generate echoes

clear
clc
close all
rng(0)

% set folder paths
config.workDir.input = '../data/input_FUS/';
config.workDir.output = '../data/output_FUS/';

% create output dir
if ~exist(config.workDir.output,'dir')
    mkdir(config.workDir.output);
    mkdir([config.workDir.output, 'Data1/']);
    mkdir([config.workDir.output, 'Data2/']);
end

% set input filenames
config.phantom.structure = 'phantomFUS.png';
config.phantom.roiParam = 'ROI_paramFUS.csv';
config.acquisition.sensitivities = 'sensitivitiesFUS.mat';

% simulation options
precontrast = 1;       % generate precontrast scans
dynamic = 1;           % generate dynamic scans
reco = 0;              % reconstruct images and save data to perflab format

%% Set input parameters for MakePhantom
config.phantom.resizeFactor = 0.8;
config.phantom.slices = 1;           % number of slices for 3D (stack of stars)

%% Set input parameters for acquisition
config.acquisition.model = 'FLASH2D';

config.acquisition.TR = 0.008;      % TR [s] (time interval between radial acquisition)
config.acquisition.TE = 0.0014447;  % TEs [s] (in case of a vector, a multi-gradient-echo acquisition assumed)
config.acquisition.FA = 15*pi/180;  % flip angle of excitation RF pulses [rad]
config.acquisition.delay = 30;      % bolus arrival time [s]
config.acquisition.krho = 0.01;     % multiplication factor accounting for spin density, gains,... here set arbitrarily, don't change (it would change the SNR which is now tuned to agree with real data)
config.acquisition.SD = 0.022;      % standard deviation of noise (0.022 set experimentally as realistic for radial sampling)
config.acquisition.r1  = 3.2;       % r1 relaxivity of tissue and blood [l/mmol/s]

%% Time axis parameters
config.acquisition.timeAxis.scanTime = 256; % scan time [s]

config.acquisition.timeAxis.Ts = config.acquisition.TR * config.phantom.slices;  % sampling period [s]
config.acquisition.timeAxis.N = ceil(config.acquisition.timeAxis.scanTime / config.acquisition.timeAxis.Ts); % number of time axis samples
Nmin = config.acquisition.timeAxis.N * 2 - 1; % minimum length for zeropadding
config.acquisition.timeAxis.NN = (2^ceil(log2(Nmin))); % the closest power of two (for fast FT)

%% Set input parameters for MakeKspaceTraj
% k-space sampling options:
%   'cartesian'
%   'radial'
%   'rosettes'
config.acquisition.kSampling.method = 'radial';

switch config.acquisition.kSampling.method
    case 'cartesian'
        % parameters for cartesian sampling
        config.acquisition.kSampling.phEncSteps = 128;
        config.acquisition.kSampling.nSamples = 128;
        config.acquisition.kSampling.repetitions = ceil(config.acquisition.timeAxis.N / config.acquisition.kSampling.phEncSteps);
    case 'radial'
        % parameters for radial sampling (using NUFFT)
        config.acquisition.kSampling.angleIncrement = 2*pi/(1+sqrt(5));    % golden angle now
        config.acquisition.kSampling.nSamples = 128;
        config.acquisition.kSampling.projections = config.acquisition.timeAxis.N;   % affects acquisition length (=Projections * TR) should be at least cca 10 min (standard duration: 10 - 15 min)
    case 'rosettes'
        % parameters for rosette sampling (using NUFFT)
        config.acquisition.kSampling.Nres = 128;            % [pixels] intended rect. matrix size of the reconstructed image
        config.acquisition.kSampling.rosettes = config.acquisition.timeAxis.N; % number of rosettes (= number of excitations)
        config.acquisition.kSampling.angleIncrement = 2*pi/(1+sqrt(5));    % golden angle now
        config.acquisition.kSampling.nSamples = 5024;       % number of samples within the rosette
        config.acquisition.rBW = 625e3;                     % [Hz] sampling frequency of acquisition
        config.acquisition.tacq0 = 0.00159;                 % [s] delay before acquisition of the 1st sample (1st TE)
        config.acquisition.kSampling.n1 = 3.5;              % number of fast oscillations
        config.acquisition.kSampling.n2 = 1.5;              % number of slow oscillations
        Tacq = config.acquisition.kSampling.nSamples / config.acquisition.rBW;       % [s] acquisition duration 
        config.acquisition.kSampling.w1 = 2*pi*config.acquisition.kSampling.n1/Tacq; % [rad/s]
        config.acquisition.kSampling.w2 = 2*pi*config.acquisition.kSampling.n2/Tacq; % [rad/s]
end

%% Set input parameters for GenerateAIF 
config.phantom.AIF.model = 'AIF_triexpG';

config.phantom.AIF.A    = 2.254; 
config.phantom.AIF.B    = 0.8053;
config.phantom.AIF.C    = 0.5381;
config.phantom.AIF.tau1 = 1.433;
config.phantom.AIF.tau2 = 2.6349;
config.phantom.AIF.tau3 = 0.07;
config.phantom.AIF.beta = 0.0847;

%% Set precontrast scan parameters
if precontrast == 1  
    model = 'IRLL';    % 'vTR' / 'vFA' / 'IRLL'
    switch model
        case 'vTR'
            pars = [15 30 50 250 500];         % TR values [ms] - conversion to [s] later
        case 'vFA'
            pars = [3 5 10 15 20 25 30];       % FA values [deg] - conversion to [rad] later
        case 'IRLL'
            config.acquisition.tau = 0.0082;   % time between excitations -> trwithoutIR       
            config.acquisition.td = 0.0121;    % time between IR and first excitation -> (TR - trwithoutIR)
            config.acquisition.tr = 0.0082;    % relaxation time at the end of excitation train -> +- tau
            config.acquisition.nProjPerInv = 1500;
            FA = 3 * pi/180; % [rad]
            TE = 0.0167621; % [s]
            pars = config.acquisition.tau;
    end
    
    % cartesian x radial x rosettes
    switch config.acquisition.kSampling.method
        case 'cartesian'
            repetitions = 20; % precontrast repetitions
        case 'radial'
            projections = 12000;
        case 'rosettes'
            rosettes = 20;
    end  
end

%% RUN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set paths
mainpath = mfilename('fullpath');
mainpath = fileparts(mainpath);
addpath(fullfile(mainpath,'PKmodels'))
addpath(fullfile(mainpath,'PKmodels/AuxFunctions'))
addpath(fullfile(mainpath,'ACQmodels'))
addpath(fullfile(mainpath,'Utils'))

% Save main dir
mainDir  = config.workDir.output;
config.outFile = [];

%% Simulate dynamic data
if dynamic == 1
    % Change main folder
    config.workDir.output = [mainDir, 'Data2/'];
    
    % Run main functions
    config = MakePhantom(config);
    config = GenerateAIF(config);
    config = PKmodeling(config); 
    config = MakeKspaceTraj(config); 
    config = GenerateEchoes(config);
    
end

%% Simulate precontrast scans
if precontrast == 1    
    config.acquisition.calibration.model = model;   % calibration type
    config.acquisition.calibration.pars = pars;     % variable parameters
    config.phantom.AIF.model = 'AIF_noCA';          % zero AIF
    
    switch config.acquisition.kSampling.method  % set repetitions/projections
        case 'cartesian'
            config.acquisition.kSampling.repetitions = repetitions;
            config.acquisition.timeAxis.N = repetitions * config.acquisition.kSampling.phEncSteps;
        case 'radial'
            config.acquisition.kSampling.projections = projections;
            config.acquisition.timeAxis.N = projections;
        case 'rosettes'
            config.acquisition.kSampling.rosettes = rosettes;
            config.acquisition.timeAxis.N = rosettes;
    end
    Nmin = config.acquisition.timeAxis.N * 2 - 1; % minimum length for zeropadding
    config.acquisition.timeAxis.NN = (2^ceil(log2(Nmin))); % the closest power of two (for fast FT)
    
    % Change main folder
    config.workDir.output = [mainDir, 'Data1/'];
    
    % Run main functions that are the same for precontrast scans
    config = MakePhantom(config);
    config = GenerateAIF(config);
    config = MakeKspaceTraj(config); 
    
    % Go through variable parameters
    for i = 1:length(pars)
        fileExt = [config.acquisition.kSampling.method, '_', ...
                    config.acquisition.calibration.model, '_', ...
                    num2str(pars(i))];
        switch model
            case 'vTR'
                config.acquisition.TR = pars(i) * 1000; % conversion to [s]
                config.outFile.table = ['tissue_', fileExt, '.mat'];
                config.outFile.echoes = ['SyntheticEchoes_', fileExt, '.mat'];
            case 'vFA'
                config.acquisition.FA = deg2rad(pars(i)); % conversion to [rad]
                config.outFile.table = ['tissue_', fileExt, '.mat'];
                config.outFile.echoes = ['SyntheticEchoes_', fileExt, '.mat'];
            case 'IRLL'
                config.acquisition.model = 'IRLL';
                config.acquisition.FA = FA;
                config.acquisition.TE = TE;
                config.outFile.table = ['tissue_', fileExt, '.mat'];
                config.outFile.echoes = ['SyntheticEchoes_', fileExt, '.mat'];
        end
        
        % Run main functions
        config = PKmodeling(config);
        config = GenerateEchoes(config);
    end
    
    % Delete output file from MakePhantom
    delete([config.workDir.output, 'tissue.mat']); 
end

%% Convert simulated data to PerfLab format
if reco
    switch config.acquisition.kSampling.method
        case 'cartesian'
            ToPerfLab(mainDir, precontrast);
        case 'radial'
            ToPerfLab_iNUFT(mainDir, precontrast, 89); % number of radials
        otherwise
            disp('No reconstruction method for selected k-sampling.')
    end
end

