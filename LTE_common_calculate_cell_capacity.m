function [...
    sector_capacity,...
    max_SINR_dB_all,...
    sector_assignment,...
    maxSINR_assignment,...
    diff_SINR_dB_all,...
    cell_sizes,...
    cell_centers...
    ] = LTE_common_calculate_cell_capacity(LTE_config,networkPathlossMap,sites,eNodeBs,varargin)
% Calculates average cell capacity based on the system settings and returns
% the cdf of the SINR for the target sector and cell (the same, just that smoother)
% (c) Josep Colom Ikuno, INTHFT, 2008
% www.nt.tuwien.ac.at

if ~isempty(varargin)
    shadow_fading_used = true;
    networkShadowFadingMap = varargin{1};
else
    shadow_fading_used = false;
end

if shadow_fading_used && LTE_config.debug_level>=1
    fprintf('Calculating average sector capacity (macroscopic and shadow fading)\n');
else
    fprintf('Calculating average sector capacity (macroscopic fading)\n');
end

%% Preallocate for the SINR matrices
N_eNodeBs = length(eNodeBs);
N_maps    = size(networkPathlossMap.pathloss,3);

% For hexagonal grid:
macro_sites     = sites(strcmp({sites.site_type},'macro'));
N_macro_eNodeBs = length([macro_sites.sectors]);

% Framework for more common macro eNodeB distribution
% macro_sector_indices = [];
% if strcmp(LTE_config.femtocells_config.mode, 'CSG')
%     for b_ = 1:num_eNodeBs % Maybe could be implemented more efficiently with cell array of 'site_type'
%         if (strcmp(eNodeBs(b_).site_type, 'macro'))
%             macro_sector_indices = [macro_sector_indices eNodeBs(b_).sectors.eNodeB_id];
%         end
%     end
% end

RRHs = [eNodeBs.RRHs];
if ~isempty(RRHs)
    RRH_parents            = [RRHs.parent_eNodeB];
    RRHs_eNodeB_tx_W       = [RRH_parents.max_power];
    RRH_parent_eNodeB_ids  = [RRH_parents.eNodeB_id];
else
    RRHs_eNodeB_tx_W = [];
    RRH_parent_eNodeB_ids          = [];
end
attached_sites         = [eNodeBs.parent_eNodeB];
site_ids_shadow_fading = [attached_sites.id max([attached_sites.id])+(1:length(RRHs))];
map_TX_W               = [eNodeBs.max_power RRHs_eNodeB_tx_W];
RX_W_map_idx           = [[eNodeBs.eNodeB_id] RRH_parent_eNodeB_ids];

% Check total number of TX antennas including RRHs
eNodeB_total_nTX       = zeros(1,N_eNodeBs);
power_ratio            = zeros(1,N_maps);
for b_=1:N_eNodeBs
    eNodeB_total_nTX(b_) = eNodeBs(b_).nTX;
    
    if ~isempty(eNodeBs(b_).RRHs)
        eNodeB_total_nTX(b_) = eNodeB_total_nTX(b_) + sum([eNodeBs(b_).RRHs.nTX]);
    end
    power_ratio(b_) = eNodeBs(b_).nTX / eNodeB_total_nTX(b_);
    
    for rrh_=1:length(eNodeBs(b_).RRHs)
        power_ratio(eNodeBs(b_).RRHs(rrh_).id) = eNodeBs(b_).RRHs(rrh_).nTX / eNodeB_total_nTX(b_);
    end
end

shadow_fading_new_map_size = [size(networkPathlossMap.pathloss,1) size(networkPathlossMap.pathloss,2)];
no_need_to_calculate_shadow_fading = ~shadow_fading_used || (shadow_fading_used && isa(networkShadowFadingMap,'channel_gain_wrappers.shadowFadingDummyMap'));
if no_need_to_calculate_shadow_fading
    shadow_fading_current_eNodeB_lin = ones(shadow_fading_new_map_size);
end

% Calculate RX power taking into account the RRHs
RX_powers_W  = zeros([size(networkPathlossMap.pathloss,1) size(networkPathlossMap.pathloss,2) length(attached_sites)]);
for m_ = 1:N_maps
    % Get shadow fading (if necessary)
    if ~no_need_to_calculate_shadow_fading
        % Resize in the log domain for precision, also avoid negative
        % values due to the interpolation
        shadow_fading_current_eNodeB_lin = max(10.^(imresize(10*log10(networkShadowFadingMap.pathloss(:,:,site_ids_shadow_fading(m_))),shadow_fading_new_map_size)/10),0);
    end
    
    TX_power_W = map_TX_W(m_);
    map_idx    = RX_W_map_idx(m_);
    
    % If RRHs present, RX power of sector-eNodeB and RRHs is combined
    % coherently
    RX_powers_W(:,:,map_idx) = RX_powers_W(:,:,map_idx) + ...
        power_ratio(m_)*TX_power_W ./ ...
        networkPathlossMap.pathloss(:,:,m_) ./ ...
        shadow_fading_current_eNodeB_lin;
end

%% Calculate SINR map for all sectors
SINR_linear_all = zeros(size(RX_powers_W));
SNR_linear_all  = zeros(size(RX_powers_W));
thermal_noise_W = 10^(LTE_config.UE.thermal_noise_density/10) / 1000 * LTE_config.bandwidth * 10^(LTE_config.UE.receiver_noise_figure/10);

tot_eNodeBs = size(RX_powers_W,3);
for s_=1:tot_eNodeBs
    SNR_linear_all(:,:,s_)  = RX_powers_W(:,:,s_) ./ (thermal_noise_W);
    SINR_linear_all(:,:,s_) = RX_powers_W(:,:,s_) ./ (sum(RX_powers_W,3) + thermal_noise_W - RX_powers_W(:,:,s_));
end
%SNR_dB_all  = 10*log10(SNR_linear_all);
SINR_dB_all = 10*log10(SINR_linear_all);
% Calculate the matrix needed to show the SINR difference map
%[SNR_dB_all_sorted,   SNR_IX] = sort(SNR_dB_all,3);
[SINR_dB_all_sorted, SINR_IX] = sort(SINR_dB_all,3);

max_SINR_dB_all    = SINR_dB_all_sorted(:,:,end);
if tot_eNodeBs>1
    diff_SINR_dB_all   = SINR_dB_all_sorted(:,:,end)-SINR_dB_all_sorted(:,:,end-1);
else
    diff_SINR_dB_all = nan(size(SINR_dB_all_sorted(:,:,1)));
end

maxSINR_assignment = SINR_IX(:,:,end);
% OSG mode: Take into account all eNodeBs.
% CSG mode: Sector assignment of femtos is superimposed a-posteriori.
if LTE_config.add_femtocells
    switch LTE_config.femtocells_config.mode
        case 'OSG'
            sector_assignment  = maxSINR_assignment;
        case 'CSG'
            SINR_dB_macro_only = SINR_dB_all(:,:,1:N_macro_eNodeBs);
            [~, SINR_IX] = sort(SINR_dB_macro_only,3);
            sector_assignment  = SINR_IX(:,:,end);
            if tot_eNodeBs > N_macro_eNodeBs % there are femto sectors
                % For each base stations determine points with r < R_indoorArea
                indoor_areas_ = (networkPathlossMap.distances < LTE_config.femtocells_config.macroscopic_pathloss_model_settings.indoorAreaRadius);
                % Get only the femto sectors:
                indoor_areas  = indoor_areas_(:,:,N_macro_eNodeBs+1:end);
                % Determine sector assignment
                [~, Femto_IX] = sort(indoor_areas,3);
                % Get indices of indoor areas and filter the rest of the map
                crop_indoor_areas = (sum(indoor_areas,3)>0);
                Femto_IX = (Femto_IX(:,:,end)+N_macro_eNodeBs).*crop_indoor_areas;
                sector_assignment = sector_assignment.*(~crop_indoor_areas) + Femto_IX;
            end
        otherwise
            error('Unknown femtocell subscriber mode');
    end
else
    sector_assignment  = maxSINR_assignment;
end

% Calculate sector sizes
cell_sizes = zeros(1,length(eNodeBs));
cell_centers_pixel = zeros(length(eNodeBs),2);
for s_idx = 1:length(eNodeBs)
    cell_sizes(s_idx) = sum(sector_assignment(:)==s_idx);
    [row,col] = find(sector_assignment==s_idx);
    cell_centers_pixel(s_idx,:) = [mean(col) mean(row)];
end

cell_centers = LTE_common_pixel_to_pos(cell_centers_pixel,networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);

SINR_dB_ROI  = max_SINR_dB_all(:);
SINR_lin_ROI = 10.^(SINR_dB_ROI/10);

%% Calculate average capacity (finally!!!)
bandwidth       = LTE_config.N_RB*LTE_config.RB_bandwidth;
CP_length_s     = LTE_config.CP_length_samples/LTE_config.fs;
symbol_length_s = LTE_config.TTI_length/(LTE_config.N_sym*2);
CP_ratio        = 1-(CP_length_s/symbol_length_s);

nTXantennas        = 1;
subcarriers_per_RB = 12;
switch nTXantennas
    case 1
        nRef_sym = 4;
    case 2
        nRef_sym = 8;
    case 4
        nRef_sym = 12;
end
subframe_size_Sym = LTE_config.N_sym*subcarriers_per_RB*2*LTE_config.N_RB;       % 2 for 2 slots (2x0.5 ms)

RefSym_ratio  = 1-(nRef_sym / (LTE_config.N_sym*subcarriers_per_RB*nTXantennas)); % Ratio of reference_symbols/total_subframe_symbols
SyncSym_ratio = 1-(72 / (subframe_size_Sym*5)); % 72 symbols used for sync every 5 subframes

% Integrate over all of the ROI (sum). Apply correction factors for used bandwidth, Cyclic Prefix and reference/sync symbols.
sector_capacity_vec     = bandwidth*CP_ratio*RefSym_ratio*SyncSym_ratio*log2(1+SINR_lin_ROI);

sector_avg_capacity_mbps = mean(sector_capacity_vec) / 1e6;
sector_min_capacity_mbps = min(sector_capacity_vec) / 1e6;
sector_max_capacity_mbps = max(sector_capacity_vec) / 1e6;

sector_capacity.avg_mbps = sector_avg_capacity_mbps;
sector_capacity.min_mbps = sector_min_capacity_mbps;
sector_capacity.max_mbps = sector_max_capacity_mbps;

end
