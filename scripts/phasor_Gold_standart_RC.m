%% Prepare meta data
clear;
close all;

% ADD DEPENDENCIES
addpath('C:\WORK_DIRECTORY\dzhigd\matlab_development\phasor_oop\classDef');
addpath('C:\WORK_DIRECTORY\dzhigd\matlab_development\phasor_oop\functions');
addpath('reconstructionPatterns');

input_param.data_path         = 'C:\WORK_DIRECTORY\dzhigd\CurrentProject\FIB_STO_BCDI\APS\data\Sample3\Sample3__165\Sample3__165_diff.mat';
% input_param.data_path         = 'C:\WORK_DIRECTORY\dzhigd\matlab_development\PHASOR_GUI_O\data_test.mat';
input_param.save_path         = 'C:\WORK_DIRECTORY\dzhigd\CurrentProject\FIB_STO_BCDI\APS\data\Sample3\Sample3__165';
input_param.mask_path         = [];
input_param.beamtime_id       = 'Dzhigaev1118';
input_param.sample_name       = 'Gold_standart_RC';
input_param.gpu               = 1;

%% Experimental parameters
input_param.energy            = 9000;   % [eV]
input_param.sample_detector_d = 0.5;   % [m]
input_param.detector_pitch    = 55e-6; % [m]
input_param.angular_step      = -0.01; % [degree]
input_param.energy_step       = [];     % [eV]

% PHASOR parameters
input_param.binning           = []; % in the detector plane
input_param.metric.type       = 'both';

input_param.partial_coherence = 0;
input_param.pc_iterations     = 20;
input_param.partial_coherence_start = 30;
input_param.partial_coherence_end = 300;

n = 1;
input_param.dataTime = getTimeStamp();

for jj = 1:n
    jj    
    % Prepare PHASOR object

    % The object has to be created
    phasor = Phasor(input_param);

    % just loading the mat file with cleaning from nans and recording data gaps as regions with -1 values
    phasor.load_mat; 

%     phasor.center_data;
    
    % create support 
    % default: autocorrelation with treshold of 0.1
    % <path> : load support from outside
    phasor.create_support();

    % create object 
    % default: constant amplitude of 1 multiplied by support
    % <type> : % amp  - random amplitude, constant phase
               % ph   - random phase, constant amplitude
               % amph - random phase, random amplitude
               % flat - constant amplitude
    phasor.create_object('ph');

    % The package can be used in various templates of scripts going below
    % %% Simple Phase retrieval pattern

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GPU part - common for any pattern %%%%%%%%%%%%%%%%%%
    if input_param.gpu
        phasor.ram2gpu;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %>>>>>>>>>>>>>>>>>>> Execute phasing pattern >>>>>>>>>>>>>>>>>>>>>>>>%
%     tripleLoop_STO_78_94;
    jesseClarckStandart
%     tripleLoop;
%     corneliusStandart;
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GPU part - common for any pattern %%%%%%%%%%%%%%%%%%
    if input_param.gpu
        phasor.gpu2ram;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Evaluation of results
    phasor.plot_metric;
    
    phasor.center_object;
    
    % Remove phase only before saving, no further phasing is possible
    phasor.remove_phase_ramp;
    
    phasor.slice3d_object;
    phasor.iso3d_object;
    % Save
%     phasor.save(string(jj));
    metric_reci(jj)  = (phasor.metric.val.reciprocal(end));
    metric_sharp(jj) = (phasor.metric.val.sharpness(end));
end

%% Experimental parameters
phasor.data_meta.energy            = 9000;   % [eV]
phasor.data_meta.sample_detector_d = 0.5;   % [m]
phasor.data_meta.detector_pitch    = 55e-6; % [m]
phasor.data_meta.angular_step      = -0.01; % [degree]
phasor.data_meta.energy_step       = [];     % [eV]

% Experimental geometry - very important for coordinate transformations
% First two dimensions in the data is the detector plane
% The rocking curve direction goes along 3rd dimension of the data
% In the case of 34idc there is a method to read marameters from .spec file
phasor.data_meta.beamline          = '34idc';
phasor.data_meta.delta             = 32.99125;
phasor.data_meta.gamma             = -11.3018;
phasor.data_meta.theta             = 0.1699971;
phasor.data_meta.chi               = 95;
phasor.data_meta.phi               = -4.03+3.03;
phasor.data_meta.mu                = 0.0;
phasor.data_meta.rocking_motor     = 'dtheta';

% phasor.transform2lab;

% figure;
% plot(metric_sharp,'-r');
% figure;
% plot(metric_reci,'-k');
% % save('phasor.data_meta.save_path\metric_reci.mat','metric_reci');