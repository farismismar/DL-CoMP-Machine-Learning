close all force;
clc;
clear all
clear global;
clear classes;

% Must be run from terminal to inherit all env variables.
% https://stackoverflow.com/questions/45733111/importing-tensorflow-in-matlab-via-python-interface

% CoMPenabled This is where the transition graph shows (plotted at the
% bottom of the file)
global CoMPDecisions;

startTime = tic;

Q = 10; % UEs per BS

% This is T_CoMP
global dComp;
dComp = 3; % TTIs.
Total_Time = 10*dComp - 1;  % 4 

global DLCoMPSINRMin;
DLCoMPSINRMin = -3.5; % dB  on average.

global epsilon;

global seed;
seed = 7;
rng(seed,'twister');

% Scenario flag.
global staticCoMP;

global model_choice;
model_choice = 'svm'; % svm or dnn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% staticCoMP = false: dynamic algorithm
%            = true: static cutoff based on DLCoMPSINRMin
staticCoMP = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (model_choice == 'svm')
    epsilon = 0.15;  % for Error
end

if (model_choice == 'dnn')
    epsilon = 0.15;  % for Error
end

simulation_type = 'tri_sector_plus_femtocells';

% Possible simulation types now:
%   - 'tri_sector'
%   - 'tri_sector_tilted', 'tri_sector_tilted_4x2', 'tri_sector_tilted_4x4'
%   - 'tri_sector_plus_femtocells'
%   - 'six_sector_tilted'
%   - 'capesso_pathlossmaps'
%   - 'omnidirectional_eNodeBs'

LTE_config = LTE_load_params(simulation_type);
%% If you want to modify something taking as a base the configuration file, do it here: here an example is show that changes the inter-eNodeB distances based on the LTE_load_params_hex_grid_tilted config file.

% Some changes to the base configuration, in case you would need/want them
LTE_config.debug_level                = 1; % basic output. % 2 is buggy.
LTE_config.bandwidth                  = 10e6; % 10 MHz
LTE_config.frequency                  = 2.1e9; % 2.1 GHz
LTE_config.channel_model.type         = 'PedA';
LTE_config.use_fast_fading            = true;
LTE_config.show_network               = 3; % show plots - Everything
LTE_config.nTX                        = 2;
LTE_config.nRX                        = 2;  % Create a 2T2R channel
LTE_config.tx_mode                    = 4;  % 4 = CLSM
LTE_config.seedRandStream             = true; % Allow reproducibility
LTE_config.RandStreamSeed             = seed;
LTE_config.scheduler                  = 'prop fair Sun'; 
LTE_config.network_source             = 'generated'; % hexagonals
LTE_config.network_geometry           = 'regular_hexagonal_grid';
LTE_config.nr_eNodeB_rings            = 1; %0; %1; % 1 ring basically means a center site plus six sites surrounding it (case 2). 0 is single macro (case 1).
LTE_config.inter_eNodeB_distance      = 100; % 100m apart.
%LTE_config.antenna_azimuth_offsett    = 0;  % Changes the reference of the azimuth at 0 degrees.
LTE_config.macroscopic_pathloss_model = 'cost231'; % Good for HN and 2100 MHz simulations.
LTE_config.macroscopic_pathloss_model_settings.environment = 'urban_macro'; %for cost231.
LTE_config.shadow_fading_type         = 'claussen';
LTE_config.shadow_fading_mean         = 0;
LTE_config.shadow_fading_sd           = 8; % 8 dB
LTE_config.eNodeB_tx_power            = 10^((46-30)/10); % 46 dBm for macro
LTE_config.site_altiude               = 0;  % average terrain height 
LTE_config.site_height                = 25; % site height above terrain
LTE_config.rx_height                  = 1.5; % UE is at 1.5 meters
LTE_config.antenna_gain_pattern       = 'kathreinTSAntenna'; % not showing as 3D antenna. 
%LTE_config.sector_azimuths            = [30 150 270];
LTE_config.antenna.electrical_downtilt= 4;
LTE_config.max_antenna_gain           = 17; % 17 dB
LTE_config.UE.thermal_noise_density   = -174; % dBm/Hz
LTE_config.cache_network              = false;
LTE_config.antenna.antenna_type = '742212';
LTE_config.antenna.frequency = 2140;
                
% Small cell layer
LTE_config.add_femtocells             = true;  % femto but configured as a pico with power
LTE_config.femtocells_config.tx_power_W = 10^((37-30)/10); % 37 dBm is 5W.
LTE_config.femtocells_config.spatial_distribution = 'homogenous density';
LTE_config.femtocells_config.femtocells_per_km2 = 100; 
%LTE_config.femtocells_config.macroscopic_pathloss_model = 'cost231'; % 'dual slope'

LTE_config.compact_results_file       = true;
LTE_config.delete_ff_trace_at_end     = true;
LTE_config.UE_cache                   = false;
LTE_config.UE.antenna_gain            = -1; % -1 dB.
LTE_config.UE.nRX                     = 2;
LTE_config.UE.receiver_noise_figure   = 7; % 7dB
%LTE_config.UE_cache_file              = 'auto';
LTE_config.adaptive_RI                = 0;
LTE_config.keep_UEs_still             = false;
LTE_config.UE_per_eNodeB              = Q;
LTE_config.UE_speed                   = 5/3.6; % Speed at which the UEs move. In meters/second: 5 Km/h = 1.38 m/s
LTE_config.map_resolution             = 20;  % 1 % the highest resolution.  WARNING: Making it too small compared to ISD will cause cells with 0 UEs!
LTE_config.pregenerated_ff_file       = 'auto';
LTE_config.trace_version              = 'v1';    % 'v1' for pregenerated precoding. 'v2' for run-time-applied precoding

LTE_config.simulation_time_tti        = Total_Time;  %LTE_config.UE.antenna_gain: should be 20*compl
%%%%%%%%%%%%%
output_results_file = LTE_sim_main(LTE_config); % This is the main line... do not re run it unless you know what you are doing.
%%%%%%%%%%%%%
simulation_data                   = load(output_results_file);
GUI_handles.aggregate_results_GUI = LTE_GUI_show_aggregate_results(simulation_data);
GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data,GUI_handles.aggregate_results_GUI);

% Generate the plot
figure(1000)
plot(0:Total_Time, [1;CoMPDecisions], 'k')
grid on;
xlabel('TTI')
ylabel('CoMP Decision')
ylim([-2,2])

if (staticCoMP)
    dec_file = 'static_comp_decisions.tikz';
else
    dec_file = sprintf('%s_comp_decisions.tikz', model_choice);
end
matlab2tikz(dec_file);

endTime = toc(startTime);
fprintf('Simulation: total time = %1.5f sec.\n', endTime);
fprintf('Simulation: quit.\n');