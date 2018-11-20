function [output_results_file, pregenerated_ff, LTE_config] = LTE_sim_main(LTE_config,varargin)
% LTE System Level simulator main simulator file.
%
% (c) Josep Colom Ikuno, INTHFT, 2008
%

% Faris
% Run it once:
% pyversion '/Library/Frameworks/Python.framework/Versions/3.6/bin/python3'

% Invoke Python code
module = py.importlib.import_module('dnn'); % a pointer to dnn.py
py.importlib.reload(module);

py.importlib.import_module('os');

global seed;
global model_choice; % are we using dnn or svm?
py.dnn.initialize_wrapper(int64(seed));
fall_back = false; % is the ML model unable to predict?  Fall back to static CoMP.
validity = false;  % is the ML model valid?
% End Faris

v = ver;
if ~isempty(strfind([v.Name],'Parallel Computing Toolbox'))
    t_ID = getCurrentTask();
    if isempty(t_ID)
        t_ID = 1;
    else
        t_ID = t_ID.ID;
    end
else
    t_ID = 1;
end


the_date = date;
if t_ID == 1
    fprintf('Vienna LTE-A System Level simulator\n');
    fprintf('(c) 2008-%s, Institute of Telecommunications (TC), TU Wien\n',the_date((end-3):end));
    fprintf(' This work has been funded by A1 Telekom Austria AG and the Christian Doppler Laboratory for Design Methodology of Signal Processing Algorithms.\n\n');
    fprintf('  By using this simulator, you agree to the license terms stated in the license agreement included with this work\n\n');
end

ticID_start_sim_begin = tic;

%% Call load_params_dependant
LTE_config  = LTE_load_params_dependant(LTE_config);
DEBUG_LEVEL = LTE_config.debug_level;
varargin{4} = []; % To avoid having an empty varargin variable

%% Load and plot BLER curves
[ BLER_curves, CQI_mapper ] = LTE_init_load_BLER_curves(LTE_config);

%% Get the eNodeBs, the Macroscopic Pathloss and the Shadow Fading
% No need to generate shadow fading maps when using network planning tool
[sites, eNodeBs, networkPathlossMap, networkShadowFadingMap LTE_config] = LTE_init_network_generation(LTE_config,varargin);
CoMP_sets = CoMP.initialize_CoMP_sets(LTE_config, eNodeBs);
varargin{2} = []; % clear varargin{2}

if networkPathlossMap.version < LTE_config.minimum_pathloss_map_version
    error('Older pathloss version found. Maybe an old cached file? If this is the case, just delete the old cache file. v%d found, v%d required',networkPathlossMap.version,LTE_config.minimum_pathloss_map_version);
end
% Compute only UEs, which are attached to an eNB within a certain inclusion range.
% eNBs determined here a-posteriori.
if isfield(LTE_config, 'compute_only_center_users')
    if LTE_config.compute_only_center_users
        LTE_config.compute_only_UEs_from_this_eNodeBs = LTE_init_determine_eNodeBs_to_compute(LTE_config, sites);
        LTE_config.default_shown_GUI_cells            = LTE_config.compute_only_UEs_from_this_eNodeBs;
    end
end

%% Plot network
if LTE_config.show_network>1
    use_subplots = true; % Change to false to have all of the plots in separate figures (useful to produce results/paper figures)
    LTE_plot_loaded_network(LTE_config,sites,eNodeBs,networkPathlossMap,CQI_mapper,use_subplots,networkShadowFadingMap);
end

%% Add a clock to each network element
networkClock = network_elements.clock(LTE_config.TTI_length);

for b_=1:length(sites)
    sites(b_).clock = networkClock;
end


if LTE_config.support_beamforming
    beamforming_coordinator = network_elements.beamforming_coordinator(LTE_config,eNodeBs);
end


%% Create users (UEs) and add schedulers
[UEs, extra_UE_cache_info, LTE_config] = LTE_init_generate_users_and_add_schedulers(LTE_config,sites,eNodeBs,networkPathlossMap,CQI_mapper,BLER_curves,networkClock);

%%
if isempty(UEs)
    warning('No UEs generated. Simulation will be skipped and no results file will be saved.');
else
    %% Generate/load the fast fading traces
    % In case of 3D channel model, a scenario-specific channel is generated in each TTI
    if ~strcmp(LTE_config.channel_model.type,'3D_UMa_fading') && ~strcmp(LTE_config.channel_model.type,'3D_UMi_fading')
        if isempty(strfind(LTE_config.pregenerated_ff_file,'.mat'))
            LTE_config.pregenerated_ff_file = sprintf('%s.mat',LTE_config.pregenerated_ff_file);
        end
        ff_file_exists = exist(LTE_config.pregenerated_ff_file,'file');
        if LTE_config.recalculate_fast_fading || (~ff_file_exists && ~LTE_config.recalculate_fast_fading)
            % Generated UE fast fading
            if DEBUG_LEVEL>=1
                fprintf('Generating UE fast fading and saving to %s\n',LTE_config.pregenerated_ff_file);
            end
            pregenerated_ff = LTE_init_get_microscale_fading_SL_trace(LTE_config);
            try
                if exist(LTE_config.pregenerated_ff_file,'file')
                    throw(MException('LTEsim:cacheExists', 'The cache file was concurrently generated during another simulation run'));
                end
                save(LTE_config.pregenerated_ff_file,'pregenerated_ff','-v7.3');
            catch err
                fprintf('Channel trace could not be saved. If needed, it will be generated again in the next run (%s).\n',err.message);
            end
        else
            % Load UE fast fading
            if LTE_config.reuse_pregenerated_ff_trace_from_last_run
                if DEBUG_LEVEL>=1
                    fprintf('Reusing UE fast fading trace from input.\n');
                end
                if isempty(varargin{1})
                    error('Extra input for LTE_sim_main needed: pregenerated_fading_trace.');
                end
                pregenerated_ff = varargin{1};
                varargin{1} = []; % clear varargin{1}
            else
                if DEBUG_LEVEL>=1
                    fprintf('Loading UE fast fading from %s \n.',LTE_config.pregenerated_ff_file);
                end
                load(LTE_config.pregenerated_ff_file,'pregenerated_ff');
                
                % v1 fading traces generate a v2 trace as master channel file,
                % then load that channel to create the v1 trace. This check
                % ensures that the originally used master file is still
                % existent, if not it throws an error. If you want to use old
                % files, or want to ignore constency, put
                % LTE_config.dont_check_consistency == 1
                % In your batch file
                if strcmp(LTE_config.trace_version, 'v1')&&...                only meaningful for v1
                        (~isfield(LTE_config, 'dont_check_consistency') ||... check if dont_check_consistency is not defined
                        (LTE_config.dont_check_consistency == 0))  %          or set to zero
                    
                    % Check if master channel file still exists. If not:
                    % Delete all v1 traces corresponding to the nTX-nRX
                    % configuration, and let the simulator calculate a fresh
                    % trace
                    if ~exist(pregenerated_ff.generated_from)
                        error(sprintf('Base channel from which this trace was generated does not exist anymore:\n %s',pregenerated_ff.generated_from))
                    end
                    
                    % Checks if master channel file has changed. If it has,
                    % delete the v1 trace and it will be recalculated based on
                    % the new channel
                    source_info = dir(pregenerated_ff.generated_from);
                    if ~isequal(pregenerated_ff.source_info , source_info)
                        error(sprintf('Base file has changed. Delete this trace to recalculate it from the new base file.\n Old base file: %s',pregenerated_ff.generated_from))
                    end
                end
            end
            
            % Wrong trace version
            switch class(pregenerated_ff)
                case 'phy_modeling.PregeneratedFastFading_v2'
                    if ~strcmp(LTE_config.trace_version,'v2')
                        error('Trace should be %s, %s found instead\n',LTE_config.trace_version,class(pregenerated_ff));
                    end
                case 'phy_modeling.PregeneratedFastFading'
                    if ~strcmp(LTE_config.trace_version,'v1')
                        error('Trace should be %s, %s found instead\n',LTE_config.trace_version,class(pregenerated_ff));
                    end
            end
            
            % Wrong number of nTX or nRX antennas
            if LTE_config.nTX~=pregenerated_ff.nTX || LTE_config.nRX~=pregenerated_ff.nRX
                error('Trace is for a %dx%d system. Config defines a %dx%d system.',pregenerated_ff.nTX,pregenerated_ff.nRX,LTE_config.nTX,LTE_config.nRX);
            end
            
            % Wrong bandwidth case
            if LTE_config.bandwidth ~= pregenerated_ff.system_bandwidth
                error('Loaded FF trace is not at the correct frequency: %3.2f MHz required, %3.2f MHz found',LTE_config.bandwidth/1e6,pregenerated_ff.system_bandwidth/1e6);
            end
            
            % Wrong UE speed case (not applicable if trace UE_speed is NaN->speed independent)
            if (pregenerated_ff.UE_speed~=LTE_config.UE_speed) && ~isnan(pregenerated_ff.UE_speed)
                error('Loaded FF trace is generated at %3.2f m/s. UE speed is %3.2f m/s. Trace cannot be used.',pregenerated_ff.UE_speed,LTE_config.UE_speed);
            end
            
            % Wrong channel type
            if ~strcmpi(pregenerated_ff.channel_type,LTE_config.channel_model.type)
                error('Loaded FF trace is not for the correct channel type. %s found, %s required.',pregenerated_ff.channel_type,LTE_config.channel_model.type);
            end
            
            % Print microscale fading trace speed
            if isnan(pregenerated_ff.UE_speed)
                if DEBUG_LEVEL>=1
                    sprintf('Microscale fading trace is speed-independent\n');
                end
            else
                if DEBUG_LEVEL>=1
                    sprintf('UE Fast fading trace at %3.2f m/s (%3.2f Km/h)\n',pregenerated_ff.UE_speed,pregenerated_ff.UE_speed*3.6);
                end
            end
            
        end
    end
    % The SINR averaging algorithm to be used. It does not make sense that
    % different UEs would use a different one, so it is the same one (saves
    % some memory)
    switch LTE_config.SINR_averaging.algorithm
        case 'MIESM'
            the_SINR_averager = utils.miesmAveragerFast(LTE_config,LTE_config.SINR_averaging.BICM_capacity_tables,LTE_config.SINR_averaging.betas);
        otherwise
            error('SINR averaging algorithm %s not supported',LTE_config.SINR_averaging.algorithm);
    end
    
    %% Some other UE initialization that requires of some of the just-calculated data
    for u_=1:length(UEs)
        
        % Add downlink channel (includes macroscopic pathloss, shadow fading and fast fading models)
        UEs(u_).downlink_channel = channel_models.downlinkChannelModel(UEs(u_));
        
        % Add eNodebs to the downlink channel (may be needed)
        UEs(u_).downlink_channel.eNodeBs = eNodeBs;
        
        % Set fast fading from the eNodeB to an attached UE.
        if strcmp(LTE_config.channel_model.type, '3D_UMa_fading') || strcmp(LTE_config.channel_model.type, '3D_UMi_fading')
            UEs(u_).downlink_channel.fast_fading_model = channel_gain_wrappers.TR36873_Fading_3D_Channel(LTE_config);
        else
            primary_DL_channel = channel_gain_wrappers.fastFadingWrapper(pregenerated_ff,'random',length(eNodeBs));
            UEs(u_).downlink_channel.fast_fading_model = primary_DL_channel;    
            
            % Set up secondary traces (only RRHs for now, and all with the same number of TX antennas)
            if LTE_config.runtime_precoding && ...
                    ~isempty(pregenerated_ff.secondary_traces) && ...
                    LTE_config.maximum_supported_secondary_channel_traces_per_UE>0
                for t_=1:LTE_config.maximum_supported_secondary_channel_traces_per_UE
                    % Change code here to set different traces to different
                    % antenna counts (i.e., traces). This would be necessary if
                    % you implement a mode where more different TX antenna
                    % counts are present
                    secondary_DL_channel = channel_gain_wrappers.fastFadingWrapper(pregenerated_ff.secondary_traces(1),'random',length(eNodeBs));
                    
                    if t_==1
                        UEs(u_).downlink_channel.secondary_fast_fading_models     = secondary_DL_channel;
                    else
                        UEs(u_).downlink_channel.secondary_fast_fading_models(t_) = secondary_DL_channel;
                    end
                end
            end
        end
        % Macroscopic pathloss
        UEs(u_).downlink_channel.set_macroscopic_pathloss_model(networkPathlossMap);
        
        % Shadow fading (data obtained from planning tools already have this information incorporated)
        if LTE_config.macroscopic_pathloss_is_model
            UEs(u_).downlink_channel.set_shadow_fading_model(networkShadowFadingMap);
        end
        
        % Cache the RB_grid object in the UE object, so as to avoid too many calls to the function. This will have to be taken into account when implementing handover
        UEs(u_).RB_grid = UEs(u_).downlink_channel.RB_grid;
        
        % Uplink channel
        UEs(u_).uplink_channel = channel_models.uplinkChannelModel(...
            UEs(u_),...
            LTE_config.N_RB,...
            LTE_config.maxStreams,...
            LTE_config.feedback_channel_delay);
        
        % Set UE SINR averaging algorithm
        UEs(u_).SINR_averager = the_SINR_averager;
    end
    
    %% Initialise the tracing
    
    % Global traces
    simulation_traces = tracing.simTraces;
    simulation_traces.TTIs_to_ignore_when_calculating_aggregates = LTE_config.feedback_channel_delay;
    % Traces from received UE feedbacks (eNodeB side)
    simulation_traces.eNodeB_rx_feedback_traces = tracing.receivedFeedbackTrace(...
        LTE_config.simulation_time_tti,...
        length(UEs),...
        LTE_config.N_RB,...
        LTE_config.maxStreams,...
        LTE_config.unquantized_CQI_feedback);
    simulation_traces.eNodeB_rx_feedback_traces.parent_results_object = simulation_traces;
    
    % eNodeB traces
    for c_=1:length(eNodeBs)
        if c_==1
            simulation_traces.eNodeB_tx_traces = tracing.enodebTrace(eNodeBs(c_).RB_grid,LTE_config.maxStreams,LTE_config.simulation_time_tti,networkClock);
        else
            simulation_traces.eNodeB_tx_traces(c_) = tracing.enodebTrace(eNodeBs(c_).RB_grid,LTE_config.maxStreams,LTE_config.simulation_time_tti,networkClock);
        end
        simulation_traces.eNodeB_tx_traces(c_).parent_results_object = simulation_traces;
        eNodeBs(c_).sector_trace   = simulation_traces.eNodeB_tx_traces(c_);
        eNodeBs(c_).feedback_trace = simulation_traces.eNodeB_rx_feedback_traces;
    end

    % UE traces
    for u_=1:length(UEs)
        UEs(u_).trace = tracing.ueTrace(...
            LTE_config.simulation_time_tti,...
            LTE_config.N_RB,...
            LTE_config.nTX,...
            LTE_config.nRX,...
            LTE_config.maxStreams,...
            LTE_config.unquantized_CQI_feedback,...
            LTE_config.trace_SINR,...
            LTE_config.scheduler_params.av_window,...
            LTE_config.TTI_length,...
            LTE_config.reduced_feedback_logs);
        UEs(u_).trace.parent_results_object = simulation_traces;
        
        if u_==1
            simulation_traces.UE_traces     = UEs(u_).trace;
        else
            simulation_traces.UE_traces(u_) = UEs(u_).trace;
        end
    end
    LTE_config.UE_count     = length(UEs);
    LTE_config.site_count   = length(sites);
    LTE_config.eNodeB_count = length(eNodeBs);
    
    %% Give the schedulers access to the UE traces
    % Then they can make decisions base on their received throughput. More
    % complex and realistic solutions may be possible, but then the eNodeB
    % should dynamically allocate resources to store UE-related data (and then
    % release once the UE is not attached to it anymore). It is easier like this :P
    for c_=1:length(eNodeBs)
        eNodeBs(c_).scheduler.set_UE_traces(simulation_traces.UE_traces);
    end
    
    %% Print all the eNodeBs
    % Print the eNodeBs (debug)
    if DEBUG_LEVEL>=2
        fprintf('eNodeB List\n');
        for b_=1:length(sites)
            sites(b_).print;
        end
        fprintf('\n');
    end
    
    %% Print all the Users
    if DEBUG_LEVEL>=2
        fprintf('User List\n');
        for u_=1:length(UEs)
            UEs(u_).print;
        end
        fprintf('\n');
    end
    
    %% Main simulation loop
    if DEBUG_LEVEL>=1 && t_ID == 1
        fprintf('Entering main simulation loop, %5.0f TTIs\n',LTE_config.simulation_time_tti);
    end
    
    % Inititialize timer
    ticID_start_sim_loop = tic;
    starting_time = toc(ticID_start_sim_loop);
    
    num_markers = 5;
    s_markings  = round(linspace(1,length(eNodeBs),num_markers));
    u_markings  = round(linspace(1,length(UEs),num_markers));
    
    % Check whether to calculate large scale fading parameters for 3D
    % channel model
    calculate_large_scale_parameters = false;
    if strcmp(LTE_config.channel_model.type, '3D_UMa_fading') || strcmp(LTE_config.channel_model.type, '3D_UMi_fading')
        calculate_large_scale_parameters = true;
    end
    
    % Faris
    global dComp;       % defined in the Main_File
    global DLCoMPSINRMin;  % defined in the Main_File
    global staticCoMP;  % defined in the Main_File
    %global useCQI;
    global CoMPDecisions;  % will collect the decisions, to be used in the final plot.
    
    data = struct; % This struct is for the whole simulation
    data.BLER = [];
    data.TBSINR = [];
    data.TBCQI = [];
    data.RX_Power_TB = [];
    
    % This struct is used for collecting measurements, eventually the
    % training data (X1, X2, ...).
    
    training_data = struct; % this is for training data
    training_data.BLER = [];
    training_data.TBSINR = [];
    training_data.TBCQI = [];
    training_data.RX_Power_TB = [];
    
    % Queues
    CoMPDecisions = [];
    % End Faris
    
    % Network clock is initialised to 0
    while networkClock.current_TTI < LTE_config.simulation_time_tti
        % First of all, advance the network clock
        networkClock.advance_1_TTI;
        if DEBUG_LEVEL>=1 && t_ID == 1
            fprintf(' TTI %5.0f/%d: ',networkClock.current_TTI,LTE_config.simulation_time_tti);
        end
            
        % Move UEs
        if ~LTE_config.keep_UEs_still
            move_all_UEs(LTE_config,UEs,networkPathlossMap,eNodeBs);
        end
        
         % After moving UEs, calculate large scale parameters for each UE
         % (autocorrelation specific for each eNodeB)
         if calculate_large_scale_parameters
             current_TTI      = networkClock.current_TTI;
             % Check whether UEs have moved
             recalculate_large_scale_parameters_for_this_eNodeBs = [];
             if current_TTI > 1
                 for jj = 1:length(UEs)
                    if UEs(jj).pos ~= simulation_traces.UE_traces(jj).position(current_TTI-1)
                       recalculate_large_scale_parameters_for_this_eNodeBs = [recalculate_large_scale_parameters_for_this_eNodeBs, UEs(jj).attached_eNodeB.eNodeB_id];
                    end
                 end
                 recalculate_large_scale_parameters_for_this_eNodeBs = unique(recalculate_large_scale_parameters_for_this_eNodeBs);
             else
                 recalculate_large_scale_parameters_for_this_eNodeBs = [eNodeBs(:).eNodeB_id];
             end
             % If UEs has moved or first TTI, calculate large scale parameters
             if ~isempty(recalculate_large_scale_parameters_for_this_eNodeBs)
                 for ii = 1:length(recalculate_large_scale_parameters_for_this_eNodeBs)
                     eNodeBs(recalculate_large_scale_parameters_for_this_eNodeBs(ii)).get_large_scale_paramters(LTE_config,networkPathlossMap);
                 end
                 if LTE_config.debug_level > 0
                     % Plot large scale parameter statistics.
                     figure_handle = figure('Name','Large scale parameter statistics');
                     utils.plotUtils.plot_large_scale_parameter_statistics(figure_handle, UEs);
                     figure_handle = figure('Name', 'Aggregate large scale parameter statistics');
                     utils.plotUtils.plot_aggregate_large_scale_parameter_statistics(figure_handle, UEs);
                 end
             end
         end
         
        if LTE_config.trace_simulation_mode
            attach_UEs_to_eNodeBs_according_to_trace(LTE_config,UEs,networkPathlossMap,eNodeBs);
        end
        
        % Beamforming
        if LTE_config.support_beamforming
           beamforming_coordinator.gather_all_channels(UEs,LTE_config.simulation_time_tti);
           [UEs,eNodeBs] = beamforming_coordinator.calculate_beamforming_vectors(UEs,eNodeBs);           
        end 
        
        % Placing the link measurement model here allows us to simulate transmissions with 0 delay
        if LTE_config.feedback_channel_delay==0 
            for u_ = 1:length(UEs)                
                % Measure SINR and prepare CQI feedback
                UEs(u_).link_quality_model(LTE_config,true);  
                if ~isempty(find(u_==u_markings,1))
                    if DEBUG_LEVEL>=1
                        fprintf('+');
                    end
                end
            end
        end
        
        % The eNodeBs receive the feedbacks from the UEs
        for s_ = 1:length(eNodeBs)
            % Receives and stores the received feedbacks from the UEs
            eNodeBs(s_).receive_UE_feedback;
        end
        
               
        % Process feedback received from UEs
        feedback_calculation.process_UE_feedback_at_eNodeB(eNodeBs, LTE_config,UEs(1).thermal_noise_W_RB);
           
        
                      
        
        % Scheduling and feedback calculation
        for s_ = 1:length(eNodeBs)
            % Scheduling
            eNodeBs(s_).schedule_users;            
                                 
                                 
            if ~isempty(find(s_==s_markings,1))
                if DEBUG_LEVEL>=1
                    fprintf('*');
                end
            end
        end 
        
              
        % For the non 0-delay case, call link quality model
        if (LTE_config.feedback_channel_delay~=0 || LTE_config.runtime_precoding)
            for u_ = 1:length(UEs)
                % Measure SINR and prepare CQI feedback
                UEs(u_).link_quality_model(LTE_config,false);
                if ~isempty(find(u_==u_markings,1))
                    if DEBUG_LEVEL>=1
                        fprintf('+');
                    end
                end
            end
        elseif ~LTE_config.runtime_precoding % zero FB delay case with non-runtime precoding
            for u_ = 1:length(UEs)
                UEs(u_).set_SINRs_according_to_schedule 
            end
        end
        
        
        
        
                
        % Call link performance model (evaluates whether TBs are received
        % corretly according to the information conveyed by the link quality
        % model. Additionally, send feedback (channel quality indicator +
        % ACK/NACK)
        for u_ = 1:length(UEs)
            UEs(u_).link_performance_model;
            UEs(u_).send_feedback;
            
            if ~isempty(find(u_==u_markings,1))
                if DEBUG_LEVEL>=1
                    fprintf('.');
                end
            end
        end
        
        if DEBUG_LEVEL>=1
            fprintf('\n');
        end
        
        if mod(networkClock.current_TTI,10)==0 && t_ID == 1
            elapsed_time = toc(ticID_start_sim_loop);
            time_per_iteration = elapsed_time / networkClock.current_TTI;
            estimated_time_to_finish = (LTE_config.simulation_time_tti - networkClock.current_TTI)*time_per_iteration;
            estimated_time_to_finish_h = floor(estimated_time_to_finish/3600);
            estimated_time_to_finish_m = estimated_time_to_finish/60 - estimated_time_to_finish_h*60;
            fprintf('Time to finish: %3.0f hours and %3.2f minutes\n',estimated_time_to_finish_h,estimated_time_to_finish_m);
        end
    
        % Faris
        % Is this dynamic CoMP?  If so, do some magic!
        close all % close all pesky figures
        if (staticCoMP == false)
            if (mod(networkClock.current_TTI,dComp) ~= 0)  % collect learning data between dComp intervals
                fprintf('CoMP Cluster: collecting DL CoMP features...\n')

                % Iterate across all UEs to see if they are eligible for DL CoMP
                % for this tti
                for jj = 1:length(UEs)
                    data.BLER = [data.BLER; simulation_traces.UE_traces(jj).BLER(1,:)];
                    data.RX_Power_TB = [data.RX_Power_TB; simulation_traces.UE_traces(jj).rx_power_tb];
                    data.TBSINR = [data.TBSINR; simulation_traces.UE_traces(jj).TB_SINR_dB(1,:)];
                    data.TBCQI= [data.TBCQI; simulation_traces.UE_traces(jj).TB_CQI(1,:)];
                    
                    training_data.BLER = [training_data.BLER; simulation_traces.UE_traces(jj).BLER(1,:)];
                    training_data.RX_Power_TB = [training_data.RX_Power_TB; simulation_traces.UE_traces(jj).rx_power_tb];
                    training_data.TBSINR = [training_data.TBSINR; simulation_traces.UE_traces(jj).TB_SINR_dB(1,:)];
                    training_data.TBCQI= [training_data.TBCQI; simulation_traces.UE_traces(jj).TB_CQI(1,:)];
                end
            else
                % Collection for dComp intervals done.  Now learn
                fprintf('CoMP Cluster: invalidated old ML model.\n');
                [validity, inverter, model] = estimateCoMPSettings(training_data, staticCoMP); % Faris: This vectorizes, cleans the measurement data, and trains ML.
                fprintf('CoMP Cluster: learning from collected DL CoMP features complete.\n')

                % Time to estimate new parameter based on the new measurements.
                newX = [];
                for jj = 1:length(UEs)   
                    newRsrp = 10*log10(sum(UEs(jj).rx_power_tb_in_current_tti)); % sum energy from all antennas
                   % if (useCQI)
                       % newCqi = sinrToCqi(UEs(jj).wideband_SINR);
                       % newX = [newCqi newRsrp];                    
                       % fprintf('CoMP Cluster: New UE reporting CQI of %d and RSRP of %0.1f dBm.\n', newX(1), newX(2));
                    %else                    
                    newX = [newX; UEs(jj).wideband_SINR newRsrp]; %#ok
                        %fprintf('CoMP Cluster: UE %d reporting an equivalent SINR of %0.1f dB and RSRP of %0.1f dBm.\n', jj, newX(1), newX(2));
                    %end
                end
                    % No -Inf due to dB/dBm
                    newX = newX(newX(:,1)>-Inf,:);
                    newX = newX(newX(:,2)>-Inf,:);

                    % No NaN
                    newX = newX(~isnan(newX(:,1)),:);
                    newX = newX(~isnan(newX(:,2)),:);
            end

            % Use the model for dComp - 1 more periods
            if (validity == true)
                % ML model is valid, use newX to find whether trigger
                % allowed.

              %  if (~isnan(sum(newX(:))) && ~isinf(newX(1)) && ~isinf(newX(2)))
                fall_back = false;
                if (model_choice == 'dnn')
                    save('newX.mat','newX');
                    label = py.dnn.predict_wrapper('newX.mat');  % RSRP then SINR in this order.
                    label = double(py.array.array('i',py.numpy.nditer(label))); % convert from Python to MATLAB
                    
                    % If AUC generates a value lower than 50%, invert label
                    if (inverter == true)
                        label = 1 - label;
                    end
                end
                
                if (model_choice == 'svm')
                    label = predict(model, newX);
                end

                % Finally, purge the data for a new learning
                training_data.BLER = [];
                training_data.RX_Power_TB = [];
                training_data.TBSINR = [];
                training_data.TBCQI = [];
            else
                % model is invalid.  Fall back to static
                fall_back = true;
                if (model_choice == 'dnn')
                    py.importlib.reload(module);  % purge the model.
                end
            end
            
            % Model is valid
            if (fall_back == false)
                %mean_vote = mean(label);
                median_vote = median(label);

                % Now based on the vote, keep things as is
                if (median_vote == 1)
                    % Trigger CoMP for the next TTI
                    LTE_config.CoMP_configuration = 'global'; % enable for the next TTI
                    CoMPDecisions = [CoMPDecisions;1];
                    fprintf('CoMP Cluster: DL CoMP decision is ENABLE.\n');                    
                else
                    % Switch off CoMP for the next TTI
                    LTE_config.CoMP_configuration = 'trival'; % disable for the next TTI
                    CoMPDecisions = [CoMPDecisions;0];
                    fprintf('CoMP Cluster: DL CoMP decision is DISABLE.\n');
                end
            end
        end
        
        % static CoMP here
        if (staticCoMP == true || fall_back == true)
            % Time to estimate new parameter based on the new measurements.
            newX = [];
            for jj = 1:length(UEs)   
                newRsrp = 10*log10(sum(UEs(jj).rx_power_tb_in_current_tti)); % sum energy from all antennas
               % if (useCQI)
                   % newCqi = sinrToCqi(UEs(jj).wideband_SINR);
                   % newX = [newCqi newRsrp];                    
                   % fprintf('CoMP Cluster: New UE reporting CQI of %d and RSRP of %0.1f dBm.\n', newX(1), newX(2));
                %else                    
                newX = [newX; UEs(jj).wideband_SINR newRsrp]; %#ok
                    %fprintf('CoMP Cluster: UE %d reporting an equivalent SINR of %0.1f dB and RSRP of %0.1f dBm.\n', jj, newX(1), newX(2));
                %end
                data.BLER = [data.BLER; simulation_traces.UE_traces(jj).BLER(1,:)];
                data.RX_Power_TB = [data.RX_Power_TB; simulation_traces.UE_traces(jj).rx_power_tb];
                data.TBSINR = [data.TBSINR; simulation_traces.UE_traces(jj).TB_SINR_dB(1,:)];
                data.TBCQI= [data.TBCQI; simulation_traces.UE_traces(jj).TB_CQI(1,:)];

            end

            % No -Inf due to dB/dBm
            newX = newX(newX(:,1)>-Inf,:);
            newX = newX(newX(:,2)>-Inf,:);

            % No NaN
            newX = newX(~isnan(newX(:,1)),:);
            newX = newX(~isnan(newX(:,2)),:);

            % Trigger DL CoMP if the trigger achieves min cutoff...
            % prctile(newX(:,1),[5,50,95])
            trigger = prctile(newX(:,1), 90)
            if (trigger >= DLCoMPSINRMin)
                LTE_config.CoMP_configuration = 'global'; % enable for the next TTI
                CoMPDecisions = [CoMPDecisions;1];
                fprintf('CoMP Cluster: DL CoMP decision is ENABLE.\n');  
            else
                LTE_config.CoMP_configuration = 'trival'; % disable for the next TTI
                CoMPDecisions = [CoMPDecisions;0];
                fprintf('CoMP Cluster: DL CoMP decision is DISABLE.\n');
            end
        end
    end % while TTI < sim time.
    % End Faris

    finish_time_s_full = toc(ticID_start_sim_begin);
    finish_time_m = floor(finish_time_s_full/60);
    finish_time_s = finish_time_s_full-finish_time_m*60;
    if DEBUG_LEVEL>=1
        % Faris
        %figure
        TBSINR = data.TBSINR(:);
        TBSINR = data.TBSINR(TBSINR ~= 0); % remove the -inf
        %cdfplot(TBSINR)
        %xlabel('SINR [dB]');
        %matlab2tikz('sinr_cdf.tikz');
        
        fprintf('Average SINR is %1.2f\n', mean(TBSINR(:)));
        data.TBCQI = data.TBCQI(~isnan(data.TBCQI(:)));
        fprintf('Average CQI is %1.0f\n', round(mean(data.TBCQI(:))));  
        data.RX_Power_TB = data.RX_Power_TB(~isinf(data.RX_Power_TB(:))); % get rid of infinity
        fprintf('Average RSRP is %1.2f\n', 10*log10(mean(data.RX_Power_TB(:))));
        fprintf('Average BLER is %1.2f%%\n', mean(data.BLER(:))*100);
        % End Faris
        
        fprintf('Simulation finished\n');
        fprintf(' Total elapsed time: %.0fm, %.0fs\n',finish_time_m,finish_time_s);
    end
    
    if ~isempty(extra_UE_cache_info)
        simulation_traces.extra_UE_info = extra_UE_cache_info;
    end
    
    if DEBUG_LEVEL>=1
        fprintf('Saving results to %s\n',LTE_config.results_file);
    end
    
    % Some post-processing
    simulation_traces.calculate_UE_aggregates;
    
    if strcmp(LTE_config.scheduler,'FFR')
        FFR_UE_mapping = eNodeBs(1).scheduler.UE_assignment;
    else
        FFR_UE_mapping = [];
    end
    
    % Compact (or not) results
    if LTE_config.compact_results_file
        all_UE_traces     = [simulation_traces.UE_traces];
        all_eNodeB_traces = [simulation_traces.eNodeB_tx_traces];
        extra_UE_info     = simulation_traces.extra_UE_info;
        
        networkPathlossMap.delete_everything_except_cell_assignments();
        
        if LTE_config.compact_results_file <=2
            for c_=1:length(all_eNodeB_traces)
                eNodeBs(c_).clear_non_basic_info;
                the_eNodeB_traces(c_).acknowledged_data = all_eNodeB_traces(c_).acknowledged_data;
                the_eNodeB_traces(c_).scheduled_RBs     = all_eNodeB_traces(c_).scheduled_RBs;
                the_eNodeB_traces(c_).RB_grid_size      = all_eNodeB_traces(c_).RB_grid_size;
                the_eNodeB_traces(c_).average_BLER      = all_eNodeB_traces(c_).average_BLER;
            end
        else
            the_eNodeB_traces = [];
        end
        
        for u_=1:length(all_UE_traces)
            UEs(u_).clear_non_basic_info();
            UEs_struct(u_) = UEs(u_).basic_information_in_struct();
            
            % Variables present in the compact<=2 mode
            if LTE_config.compact_results_file <= 2
                the_UE_traces(u_).RI              = all_UE_traces(u_).RI;
                the_UE_traces(u_).nCodewords      = all_UE_traces(u_).nCodewords;
                the_UE_traces(u_).position        = all_UE_traces(u_).position;
                the_UE_traces(u_).attached_site   = all_UE_traces(u_).attached_site;
                the_UE_traces(u_).attached_eNodeB = all_UE_traces(u_).attached_eNodeB;
                the_UE_traces(u_).wideband_SINR   = all_UE_traces(u_).wideband_SINR;
                the_UE_traces(u_).SNR_dB          = all_UE_traces(u_).SNR;
                the_UE_traces(u_).UE_was_disabled = all_UE_traces(u_).UE_was_disabled;
                the_UE_traces(u_).SINR            = all_UE_traces(u_).SINR;
                
            end
            
            % Variables present in the compact<=1 mode
            if LTE_config.compact_results_file <= 1
                the_UE_traces(u_).attached_site_type  = all_UE_traces(u_).attached_site_type;
                the_UE_traces(u_).TB_size             = all_UE_traces(u_).TB_size;
                the_UE_traces(u_).ACK                 = all_UE_traces(u_).ACK;
                the_UE_traces(u_).TB_CQI              = all_UE_traces(u_).TB_CQI;
                the_UE_traces(u_).CQI_sent            = all_UE_traces(u_).CQI_sent;
                the_UE_traces(u_).TB_SINR_dB          = all_UE_traces(u_).TB_SINR_dB;
				the_UE_traces(u_).RX_Power_TB         = all_UE_traces(u_).rx_power_tb;
                the_UE_traces(u_).rx_power_tb         = all_UE_traces(u_).rx_power_tb;
                the_UE_traces(u_).interference_powers = all_UE_traces(u_).interference_powers;
                the_UE_traces(u_).throughput_Mbps     = all_UE_traces(u_).throughput_Mbps;
                the_UE_traces(u_).spectral_efficiency_bit_per_cu = all_UE_traces(u_).spectral_efficiency_bit_per_cu;
                the_UE_traces(u_).channel_uses                   = all_UE_traces(u_).channel_uses;   
            end
            
            the_UE_traces(u_).average_throughput_Mbps                = all_UE_traces(u_).average_throughput_Mbps;
            the_UE_traces(u_).average_spectral_efficiency_bit_per_cu = all_UE_traces(u_).average_spectral_efficiency_bit_per_cu;
            the_UE_traces(u_).average_energy_per_bit                 = all_UE_traces(u_).average_energy_per_bit;
            the_UE_traces(u_).average_RBs_per_TTI                    = all_UE_traces(u_).average_RBs_per_TTI;
            the_UE_traces(u_).mean_wideband_SINR                     = mean(all_UE_traces(u_).wideband_SINR);           
        end
        
        UEs_bak             = UEs;
        eNodeBs_bak         = sites;
        eNodeBs_sectors_bak = eNodeBs;
        if LTE_config.compact_results_file<=1
            UEs = UEs_struct; % Saves space
        else
            UEs                = [];
            networkPathlossMap = [];
            sites       = [];
            eNodeBs            = [];
        end
        save(fullfile(LTE_config.results_folder,LTE_config.results_file),'LTE_config','networkPathlossMap','sites','eNodeBs','UEs','FFR_UE_mapping','the_UE_traces','the_eNodeB_traces','extra_UE_info','finish_time_s_full','-v7.3');
 
        UEs             = UEs_bak;
        sites         = eNodeBs_bak;
        eNodeBs = eNodeBs_sectors_bak;
        
        if DEBUG_LEVEL>=1
            fprintf('Only UE some UE and eNodeB traces saved (compact results file)\n');
        end
    else
        % Some options to save space
        if LTE_config.delete_ff_trace_at_end
            pregenerated_ff = [];
        end
        
        if LTE_config.delete_pathloss_at_end
            networkPathlossMap.pathloss = [];
            if exist('networkShadowFadingMap','var')
                networkPathlossMap.pathloss = [];
            else
                % Do nothing
            end
        end
        save(fullfile(LTE_config.results_folder,LTE_config.results_file),'-v7.3');
        if DEBUG_LEVEL>=1
            fprintf('Traces saved in the standard format\n');
        end
    end
end

output_results_file = fullfile(LTE_config.results_folder,LTE_config.results_file);

% Help Matlab clean the memory
utils.miscUtils.tidy_up_memory_before_closing(UEs,eNodeBs,sites);

function attach_UEs_to_eNodeBs_according_to_trace(LTE_config,UEs,networkPathlossMap,eNodeBs)
% Attaches each UE to a cell according to the data in the trace
if ~isempty(UEs)
    current_TTI = UEs(1).clock.current_TTI;
    
    for u_ = 1:length(UEs)
        UE_trace               = UEs(u_).downlink_channel.macroscopic_pathloss_model.pathloss(u_);
        current_trace_time     = floor((current_TTI-1)/LTE_config.TTI_per_trace_step)+1;
        current_trace_position = find(UE_trace.time_idxs==current_trace_time,1,'first');
        UE_is_active           = ~isempty(current_trace_position) && networkPathlossMap.pathloss(u_).attached_cell(current_trace_position)~=0;
        
        UEs(u_).deactivate_UE = ~UE_is_active;
        
        if isempty(UEs(u_).attached_eNodeB) && UE_is_active
            % Attach the UE
            eNodeB_to_be_attached = networkPathlossMap.pathloss(u_).attached_cell(current_trace_position);
            if eNodeB_to_be_attached~=0
                eNodeBs(eNodeB_to_be_attached).attachUser(UEs(u_));
                
                % Set a new channel realization just in case we jump from unattached to attached state
                pregenerated_ff = UEs(u_).downlink_channel.fast_fading_model.ff_trace;
                N_eNodeBs       = length(UEs(u_).downlink_channel.fast_fading_model.interfering_starting_points);
                UEs(u_).downlink_channel.set_fast_fading_model_model(channel_gain_wrappers.fastFadingWrapper(pregenerated_ff,'random',N_eNodeBs));
            end
        elseif ~isempty(UEs(u_).attached_eNodeB) && UE_is_active
            % Check if we need to change it
            eNodeB_to_be_attached = networkPathlossMap.pathloss(u_).attached_cell(current_trace_position);
            if eNodeB_to_be_attached~=0 && UEs(u_).attached_eNodeB.eNodeB_id~=eNodeB_to_be_attached
                UEs(u_).start_handover(eNodeBs(eNodeB_to_be_attached));
            end
        elseif ~isempty(UEs(u_).attached_eNodeB) && ~UE_is_active
            % Deattach the UE
            UEs(u_).attached_eNodeB.deattachUser(UEs(u_));
        end
    end
end

function move_all_UEs(LTE_config,UEs,networkPathlossMap,eNodeBs)
% This function moves each UE and for such cases in which a given UE is moved
% outside of the ROI, it places it in a random point in the ROI.

some_UE_out_of_ROI_this_TTI = false;
for u_ = 1:length(UEs)
    
       
    % If a UE went outside of ROI, relocate him somewhere else. Beam me up, Scotty! Take me somewhere in the map!!
    %
    % It also now implements a very simple handover procedure meant more as
    % an example so as to how to implement such a feature than a real handover.
    
    some_UE_out_of_ROI_this_TTI = false;
        
    %save old UE position
    old_pos     = UEs(u_).pos;
    old_pos_pix = LTE_common_pos_to_pixel([old_pos(1) old_pos(2)],networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);
    old_eNodeB_id  = networkPathlossMap.sector_assignment(old_pos_pix(2),old_pos_pix(1));
    
    %move UE
    UEs(u_).move;
    [x_range, y_range] = networkPathlossMap.valid_range;
    ROI_teleport       = ~UEs(u_).is_in_roi(x_range,y_range);
    
    %save new UE position
    new_pos     = UEs(u_).pos;
    new_pos_pix = LTE_common_pos_to_pixel([new_pos(1) new_pos(2)],networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);
    if ~ROI_teleport
        new_eNodeB_id  = networkPathlossMap.sector_assignment(new_pos_pix(2),new_pos_pix(1));
    end
    
    
    if ~LTE_config.support_handover        
        
        %check whether sector has changed
        if new_eNodeB_id ~= old_eNodeB_id || ROI_teleport
            ind = zeros(length(find(networkPathlossMap.sector_assignment == old_eNodeB_id)),2);
            [ind(:,1) ind(:,2)] = find(networkPathlossMap.sector_assignment == old_eNodeB_id);
            new_pix             = randi(length(ind(:,1)),1);
            UEs(u_).pos         = LTE_common_pixel_to_pos([ind(new_pix,2) ind(new_pix,1)],networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);
        end       
    else     
        
%         handover_requested = UEs(u_).cell_change.requested;
        if new_eNodeB_id ~= old_eNodeB_id || ROI_teleport
            if ROI_teleport
                new_pos      = networkPathlossMap.random_position;
                UEs(u_).pos  = new_pos;
            
                % Actually it should not be done like this. Measure all the neighboring cells' SINR and then decide which one is better
                [new_site_id, new_sector_id, new_eNodeB_id] = networkPathlossMap.cell_assignment(new_pos);
            end
        
            % Deattach UE from old eNodeB and reattach to new one
            UEs(u_).start_handover(eNodeBs(new_eNodeB_id));
        
            
            
            % Accordingly activate or deactivate UE
            % Check whether this UE should be deativated to speed-up
            % simulation. Allows for run-time changes.
            if ~isempty(LTE_config.compute_only_UEs_from_this_eNodeBs)
                if isempty(find(UEs(u_).attached_eNodeB.eNodeB_id==LTE_config.compute_only_UEs_from_this_eNodeBs,1))
                    % Deactivate UE
                    UEs(u_).deactivate_UE = true;
                else
                    % Activate UE (already activated by default, but just in case)
                    UEs(u_).deactivate_UE = false;
                end
            end
        
            
            % Print some debug
            if ~some_UE_out_of_ROI_this_TTI
                if LTE_config.debug_level>=1
                    fprintf(1,'\n');
                end
                some_UE_out_of_ROI_this_TTI = true;
            end
            if LTE_config.debug_level>=1
                if ROI_teleport
                    fprintf('              UE %g going out of ROI, teleporting to %g %g. eNodeB %g -> eNodeB %g\n',UEs(u_).id,new_pos(1),new_pos(2),old_eNodeB_id,new_eNodeB_id);
                else
                    fprintf('              UE %g handover request. eNodeB %g -> eNodeB %g\n',UEs(u_).id,old_eNodeB_id,new_eNodeB_id);
                end
            end
        end
    end
end
if some_UE_out_of_ROI_this_TTI
    if LTE_config.debug_level>=1
        fprintf('              ');
    end
end
