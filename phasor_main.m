%% Prepare meta data
clear;
close all;

addpath('classDef');
addpath('functions');
addpath('reconstructionPatterns');

input_param.data_path         = 'E:\DDzhigaev\WorkFolder\Projects\Flexo-photovoltaics\Experiments\APS_34IDC_2018_run3\analysis\processed\Sample3__78  82  86  90  94\data.mat';
% input_param.data_path         = 'data_test.mat';
input_param.save_folder       = fullfile('/data','netapp','dzhigd','Experiments','34_IDC_APS_STO_AFM_run3_2018','analysis','processed');
input_param.white_field_path  = fullfile('/data/netapp/dzhigd/Experiments/34_IDC_APS_STO_AFM_run3_2018/APS_34IDC_viewer','white_field_1.mat');
input_param.dark_field_path   = fullfile('/data/netapp/dzhigd/Experiments/34_IDC_APS_STO_AFM_run3_2018/APS_34IDC_viewer','dark_field_1.mat');
input_param.mask_path         = '/data/netapp/dzhigd/Experiments/34_IDC_APS_STO_AFM_run3_2018/analysis/processed/Sample1_B_pristine_mask.mat';
input_param.beamtime_id       = 'Dzhigaev1118';
input_param.sample_name       = 'test';
input_param.gpu               = 1;

% Experimental parameters
input_param.energy            = 8840;   % [eV]
input_param.sample_detector_d = 1.86;   % [m]
input_param.detector_pitch    = 172e-6; % [m]
input_param.angular_step      = 0.0025; % [degree]
input_param.energy_step       = [];     % [eV]

% Experimental geometry - very important for coordinate transformations
% init_param

% PHASOR parameters
input_param.binning           = 1;
input_param.metric.type       = 'reciprocal';
input_param.partial_coherence = 1;

%% Prepare PHASOR object

% The object has to be created
ph = Phasor(input_param);

% just loading the mat file with cleaning from nans and recording data gaps as regions with -1 values
ph.load_mat; 

ph.center_data;

% create support 
% default: autocorrelation with treshold of 0.1
% <path> : load support from outside
ph.create_support();

% create object 
% default: constant amplitude of 1 multiplied by support
% <type> : % amp  - random amplitude, constant phase
           % ph   - random phase, constant amplitude
           % amph - random phase, random amplitude
           % flat - constant amplitude
ph.create_object('amph');

%% The package can be used in various templates of scripts going below
% %% Simple Phase retrieval pattern

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPU part - common for any pattern %%%%%%%%%%%%%%%%%%
if input_param.gpu
    ph.ram2gpu;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%>>>>>>>>>>>>>>>>>>> Execute phasing pattern >>>>>>>>>>>>>>>>>>>>>>>>%
tripleLoop;
% jesseClarckStandart
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPU part - common for any pattern %%%%%%%%%%%%%%%%%%
if input_param.gpu
    ph.gpu2ram;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Evaluation of results
ph.plot_metric;

% Remove phase only before saving, no further phasing is possible
ph.remove_phase_ramp;
ph.slice3d;
ph.iso3d;

%% Save
