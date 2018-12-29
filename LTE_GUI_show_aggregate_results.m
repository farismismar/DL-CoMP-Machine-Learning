function varargout = LTE_GUI_show_aggregate_results(varargin)
% LTE_GUI_SHOW_AGGREGATE_RESULTS MATLAB code for LTE_GUI_show_aggregate_results.fig
%      LTE_GUI_SHOW_AGGREGATE_RESULTS, by itself, creates a new LTE_GUI_SHOW_AGGREGATE_RESULTS or raises the existing
%      singleton*.
%
%      H = LTE_GUI_SHOW_AGGREGATE_RESULTS returns the handle to a new LTE_GUI_SHOW_AGGREGATE_RESULTS or the handle to
%      the existing singleton*.
%
%      LTE_GUI_SHOW_AGGREGATE_RESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LTE_GUI_SHOW_AGGREGATE_RESULTS.M with the given input arguments.
%
%      LTE_GUI_SHOW_AGGREGATE_RESULTS('Property','Value',...) creates a new LTE_GUI_SHOW_AGGREGATE_RESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LTE_GUI_show_aggregate_results_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LTE_GUI_show_aggregate_results_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LTE_GUI_show_aggregate_results

% Last Modified by GUIDE v2.5 04-Apr-2012 17:05:02

% Plots the cell aggregate results for the selected eNodeBs
%
% (c) Josep Colom Ikuno, INTHFT, 2012

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LTE_GUI_show_aggregate_results_OpeningFcn, ...
                   'gui_OutputFcn',  @LTE_GUI_show_aggregate_results_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before LTE_GUI_show_aggregate_results is made visible.
function LTE_GUI_show_aggregate_results_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LTE_GUI_show_aggregate_results (see VARARGIN)

% Choose default command line output for LTE_GUI_show_aggregate_results
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Save user data (the simulation results)
simulation_data = varargin{1};
set(hObject,'UserData',simulation_data);

UE_list = cell(1,length(simulation_data.UEs));
for u_=1:length(simulation_data.UEs)
    UE_list{u_} = sprintf('%g',u_);
end

cell_list = cell(1,length(simulation_data.eNodeBs));
for c_=1:length(simulation_data.eNodeBs)
    cell_list{c_} = sprintf('%g',c_);
end

simulation_data.UE_list   = UE_list;
simulation_data.cell_list = cell_list;

% Initialize checkboxes
set(handles.cell_all,'Value',true);
set(handles.scatterplot_mean_checkbox,'Value',true);
set(handles.scatterplot_points_checkbox,'Value',true);
set(handles.constant_axes_limits_checkbox,'Value',false);

% Fill list boxes
set(handles.cell_listbox,'String',cell_list);
set(handles.cell_listbox,'Max',length(simulation_data.eNodeBs));
set(handles.cell_listbox,'Min',0);

if isfield(simulation_data.LTE_config,'default_shown_GUI_cells') && ~isempty(simulation_data.LTE_config.default_shown_GUI_cells)
    cells_to_plot = simulation_data.LTE_config.default_shown_GUI_cells;
    if simulation_data.LTE_config.compact_results_file
        cells_to_plot = cells_to_plot(cells_to_plot<=length(simulation_data.the_eNodeB_traces)); % Filtero out possible out of range values
    else
        cells_to_plot = cells_to_plot(cells_to_plot<=length(simulation_data.simulation_traces.eNodeB_tx_traces)); % Filtero out possible out of range values
    end
    set(handles.cell_listbox,'Value',cells_to_plot);
else
    set(handles.cell_listbox,'Value',1:length(simulation_data.eNodeBs));
end

% Plot main plot
plot_aggregate_results(handles)

% UIWAIT makes LTE_GUI_show_aggregate_results wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = LTE_GUI_show_aggregate_results_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in cell_listbox.
function cell_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to cell_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cell_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cell_listbox
set(handles.cell_all,'Value',false);
set(handles.cell_none,'Value',false);
plot_aggregate_results(handles);


% --- Executes during object creation, after setting all properties.
function cell_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open_plots_in_new_figure.
function open_plots_in_new_figure_Callback(hObject, eventdata, handles)
% hObject    handle to open_plots_in_new_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_aggregate_results(handles,false);


% --- Executes on button press in cell_all.
function cell_all_Callback(hObject, eventdata, handles)
% hObject    handle to cell_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_all
simulation_data = get(handles.figure1,'UserData');
if get(hObject,'Value')
    set(handles.cell_none,'Value',false);
    set(handles.cell_listbox,'Value',1:length(simulation_data.eNodeBs));
    plot_aggregate_results(handles);
end


% --- Executes on button press in cell_none.
function cell_none_Callback(hObject, eventdata, handles)
% hObject    handle to cell_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_none
simulation_data = get(handles.figure1,'UserData');
if get(hObject,'Value')
    set(handles.cell_all,'Value',false);
    set(handles.cell_listbox,'Value',[]);
    plot_aggregate_results(handles);
end

function plot_aggregate_results(handles,varargin)

% Faris
global staticCoMP;
global model_choice;
% End Faris

if isempty(varargin)
    plot_in_GUI = true;
else
    plot_in_GUI = varargin{1};
end

% Load data
simulation_data        = get(handles.figure1,'UserData');
eNodeB_sites           = simulation_data.sites;
eNodeBs                = simulation_data.eNodeBs;
UEs                    = simulation_data.UEs;
networkPathlossMap     = simulation_data.networkPathlossMap;
cells_to_plot          = get(handles.cell_listbox,'Value');
constant_axes_limits   = get(handles.constant_axes_limits_checkbox,'Value');
plot_mean_scatter      = get(handles.scatterplot_mean_checkbox,'Value');
plot_points_scatter    = get(handles.scatterplot_points_checkbox,'Value');

if isfield(simulation_data,'simulation_traces')
    % non-compact results
    the_UE_traces     = [simulation_data.simulation_traces.UE_traces];
    the_eNodeB_traces = [simulation_data.simulation_traces.eNodeB_tx_traces];
else
    the_UE_traces     = simulation_data.the_UE_traces;
    the_eNodeB_traces = simulation_data.the_eNodeB_traces;
end
cells_to_plot         = cells_to_plot(cells_to_plot<=length(the_eNodeB_traces)); % Filtero out possible out of range values
N_UEs = length(UEs);

if plot_in_GUI
    UE_throughput_axes           = handles.UE_throughput_axes;
    UE_spectral_eff_axes         = handles.UE_spectral_eff_axes;
    UE_wideband_SINR_axes        = handles.UE_wideband_SINR_axes;
    UE_SINR_to_throughput_axes   = handles.UE_SINR_to_throughput_axes;
    UE_SINR_to_spectral_eff_axes = handles.UE_SINR_to_spectral_eff_axes;
else
    a_figure                     = figure;
    UE_throughput_axes           = axes('Parent',a_figure);
    a_figure                     = figure;
    UE_spectral_eff_axes         = axes('Parent',a_figure);
    a_figure                     = figure;
    UE_wideband_SINR_axes        = axes('Parent',a_figure);
    a_figure                     = figure;
    UE_SINR_to_throughput_axes   = axes('Parent',a_figure);
    a_figure                     = figure;
    UE_SINR_to_spectral_eff_axes = axes('Parent',a_figure);
end

% Clear old plots
cla(UE_throughput_axes);
cla(UE_spectral_eff_axes);
cla(UE_wideband_SINR_axes);
cla(UE_SINR_to_throughput_axes);
cla(UE_SINR_to_spectral_eff_axes);
set(handles.cell_statistics_text,'String',[]);

if ~isempty(cells_to_plot)
    [UEs_to_use cell_sum_throughput] = utils.resultsFileReader.get_UEs_in_given_cells(cells_to_plot,the_UE_traces);
    
    the_UE_traces_to_plot = the_UE_traces(UEs_to_use);
    wideband_SINRs_all = reshape([the_UE_traces.wideband_SINR],simulation_data.LTE_config.simulation_time_tti,[]);
    wideband_SINRs_all = wideband_SINRs_all(end, :);%mean(wideband_SINRs_all, 1);
    
    % To get the average RB occupancy
    nRB                     = unique([the_eNodeB_traces(cells_to_plot).RB_grid_size]); % It should be a unique value!
    RB_occupancy_ratio      = double( [the_eNodeB_traces(cells_to_plot).scheduled_RBs])/nRB;
    mean_RB_occupancy_ratio = mean(RB_occupancy_ratio);
    
    % To get the axes limits
    throughput_Mbps_ecdf_all = utils.miscUtils.ecdf([the_UE_traces.average_throughput_Mbps]);
    spectral_eff_ecdf_all    = utils.miscUtils.ecdf([the_UE_traces.average_spectral_efficiency_bit_per_cu]);
    wideband_SINR_ecdf_all   = utils.miscUtils.ecdf(wideband_SINRs_all);
    
    % The actual values
    throughput_Mbps_ecdf = utils.miscUtils.ecdf([the_UE_traces_to_plot.average_throughput_Mbps]);
    fairness_index       = sum(throughput_Mbps_ecdf.input_data).^2 / sum(throughput_Mbps_ecdf.input_data.^2) / sum(isfinite(throughput_Mbps_ecdf.input_data));
    spectral_eff_ecdf    = utils.miscUtils.ecdf([the_UE_traces_to_plot.average_spectral_efficiency_bit_per_cu]);
    wideband_SINR_ecdf   = utils.miscUtils.ecdf(wideband_SINRs_all(UEs_to_use));
    
    % Wideband SINR-to-throughput plot
    wideband_SINR_vector                       = wideband_SINR_ecdf.input_data;
    spectral_eff_vector                        = spectral_eff_ecdf.input_data;
    throughput_vector                          = throughput_Mbps_ecdf.input_data;
    [wideband_SINR_binned  throughput_binned   numel_bins] = utils.miscUtils.fit_scatterplot_data(wideband_SINR_vector,throughput_vector,50);
    [wideband_SINR_binned2 spectral_eff_binned numel_bins] = utils.miscUtils.fit_scatterplot_data(wideband_SINR_vector,spectral_eff_vector,50);
    
    % Throughput ECDF
    hold(UE_throughput_axes,'all');
    plot(UE_throughput_axes,throughput_Mbps_ecdf.x,throughput_Mbps_ecdf.f,'blue');
    grid(UE_throughput_axes,'on');
    title(UE_throughput_axes,'UE average throughput');
    xlabel(UE_throughput_axes,'average UE throughput [Mb/s]');
    ylabel(UE_throughput_axes,'F(x)');
    if constant_axes_limits
        xlim(UE_throughput_axes,[throughput_Mbps_ecdf_all.min throughput_Mbps_ecdf_all.max]);
    else
        xlim(UE_throughput_axes,'auto');
    end
    hold(UE_throughput_axes,'off');
    
    % Spectral efficiency ECDF
    hold(UE_spectral_eff_axes,'all');
    plot(UE_spectral_eff_axes,spectral_eff_ecdf.x,spectral_eff_ecdf.f,'blue');
    grid(UE_spectral_eff_axes,'on');
    title(UE_spectral_eff_axes,'UE average spectral efficiency');
    xlabel(UE_spectral_eff_axes,'average UE spectral efficiency [bit/cu]');
    ylabel(UE_spectral_eff_axes,'F(x)');
    if constant_axes_limits
        xlim(UE_spectral_eff_axes,[spectral_eff_ecdf_all.min spectral_eff_ecdf_all.max]);
    else
        xlim(UE_spectral_eff_axes,'auto');
    end
    hold(UE_spectral_eff_axes,'off');
    
    % Wideband SINR ECDF
    hold(UE_wideband_SINR_axes,'all');
    plot(UE_wideband_SINR_axes,wideband_SINR_ecdf.x,wideband_SINR_ecdf.f,'blue');   
    grid(UE_wideband_SINR_axes,'on');
    title(UE_wideband_SINR_axes,'UE wideband SINR');
    xlabel(UE_wideband_SINR_axes,'UE wideband SINR [dB]');
    ylabel(UE_wideband_SINR_axes,'F(x)');
    if constant_axes_limits
        xlim(UE_wideband_SINR_axes,[wideband_SINR_ecdf_all.min wideband_SINR_ecdf_all.max]);
    else
        xlim(UE_wideband_SINR_axes,'auto');
    end
    hold(UE_wideband_SINR_axes,'off');
    
    % SINR-to-throughput
    hold(UE_SINR_to_throughput_axes,'all');
    if plot_points_scatter
        if plot_in_GUI
            scatter(UE_SINR_to_throughput_axes,wideband_SINR_vector,throughput_vector,'.b');
        else
            if plot_mean_scatter
                scatter(UE_SINR_to_throughput_axes,wideband_SINR_vector,throughput_vector,'.b','SizeData',50);
            else
                scatter(UE_SINR_to_throughput_axes,wideband_SINR_vector,throughput_vector,'.b','SizeData',150);
            end
        end
    end
    if plot_mean_scatter
        if plot_in_GUI
            scatter(UE_SINR_to_throughput_axes,wideband_SINR_binned,throughput_binned,'.r');
        else
            scatter(UE_SINR_to_throughput_axes,wideband_SINR_binned,throughput_binned,'.r','SizeData',350);
        end
    end
    grid(UE_SINR_to_throughput_axes,'on');
    title(UE_SINR_to_throughput_axes,'UE wideband SINR-to-throughput mapping');
    xlabel(UE_SINR_to_throughput_axes,'UE wideband SINR [dB]');
    ylabel(UE_SINR_to_throughput_axes,'average UE throughput [Mb/s]');
    if constant_axes_limits
        xlim(UE_SINR_to_throughput_axes,[wideband_SINR_ecdf_all.min wideband_SINR_ecdf_all.max]);
        ylim(UE_SINR_to_throughput_axes,[throughput_Mbps_ecdf_all.min throughput_Mbps_ecdf_all.max]);
    else
        xlim(UE_SINR_to_throughput_axes,'auto');
        ylim(UE_SINR_to_throughput_axes,'auto');
    end
    hold(UE_SINR_to_throughput_axes,'off');
    
    % Faris
    if (staticCoMP)
        dec_file = 'static_average_ue_througput_cdf.tikz';
    else
        dec_file = sprintf('%s_average_ue_througput_cdf.tikz', model_choice);
    end
    matlab2tikz(dec_file);    
    % End Faris
    
    % SINR-to-spectral efficiency
    hold(UE_SINR_to_spectral_eff_axes,'all');
    if plot_points_scatter
        if plot_in_GUI
            scatter(UE_SINR_to_spectral_eff_axes,wideband_SINR_vector,spectral_eff_vector,'.b');
        else
            if plot_mean_scatter
                scatter(UE_SINR_to_spectral_eff_axes,wideband_SINR_vector,spectral_eff_vector,'.b','SizeData',50);
            else
                scatter(UE_SINR_to_spectral_eff_axes,wideband_SINR_vector,spectral_eff_vector,'.b','SizeData',150);
            end
        end
    end
    if plot_mean_scatter
        if plot_in_GUI
            scatter(UE_SINR_to_spectral_eff_axes,wideband_SINR_binned,spectral_eff_binned,'.r');
        else
            scatter(UE_SINR_to_spectral_eff_axes,wideband_SINR_binned,spectral_eff_binned,'.r','SizeData',350);
        end
    end
    grid(UE_SINR_to_spectral_eff_axes,'on');
    title(UE_SINR_to_spectral_eff_axes,'UE wideband SINR-to-spectral efficiency mapping');
    xlabel(UE_SINR_to_spectral_eff_axes,'UE wideband SINR [dB]');
    ylabel(UE_SINR_to_spectral_eff_axes,'average UE specctral efficiency [bit/cu]');
    if constant_axes_limits
        xlim(UE_SINR_to_spectral_eff_axes,[wideband_SINR_ecdf_all.min wideband_SINR_ecdf_all.max]);
        ylim(UE_SINR_to_spectral_eff_axes,[spectral_eff_ecdf.min spectral_eff_ecdf.max]);
    else
        xlim(UE_SINR_to_spectral_eff_axes,'auto');
        ylim(UE_SINR_to_spectral_eff_axes,'auto');
    end
    hold(UE_SINR_to_spectral_eff_axes,'off');
    
    cell_sum_throughput_non_NaN      = cell_sum_throughput(isfinite(cell_sum_throughput));
    cell_sum_throughput_non_NaN_mean = mean(cell_sum_throughput_non_NaN);
    ignored_cells = sum(isnan(cell_sum_throughput));
    
    sep = '---------------------------------------';
    statistics_text = sprintf([
        '%s\nSimulations statistics:\n\n',...
        '%g cells, %g UEs\n',...
        'Simulation length: %g TTIs\n',...
        'Scheduler: %s\n',...
        'Mode: %gx%g, %s\n',...
        '%s\n',...
        'Cell statistics:\n\n',...
        'Fairness index: %g\n',...
        'Peak/Avg/Edge UE throughput:\n',...
        '%3.2f/%3.2f/%3.2f Mb/s\n',...
        'Average cell throughput: %3.2fMb/s\n',...
        'Ignored cells (disabled): %g\n',...
        'mean RB occupancy: %3.2f%%\n%s',...
        ],...
        sep,...
        length(cells_to_plot),...
        sum(UEs_to_use),...
        simulation_data.LTE_config.simulation_time_tti,...
        simulation_data.LTE_config.scheduler,...
        simulation_data.LTE_config.nTX,...
        simulation_data.LTE_config.nRX,...
        utils.miscUtils.tx_mode_to_string(simulation_data.LTE_config.tx_mode),...
        sep,...
        fairness_index,...
        throughput_Mbps_ecdf.p95,...
        throughput_Mbps_ecdf.mean_x,...
        throughput_Mbps_ecdf.p05,...
        cell_sum_throughput_non_NaN_mean,...
        ignored_cells,...
        mean_RB_occupancy_ratio*100,sep);
    if plot_in_GUI
        set(handles.cell_statistics_text,'String',statistics_text);
    else
        set(handles.cell_statistics_text,'String',statistics_text);
        fprintf('%s\n',statistics_text);
    end
    
    % Faris
    fprintf('Overall peak/average/edge throughput: %3.2f/%3.2f/%3.2f Mb/s\n.', throughput_Mbps_ecdf.p95,...
        throughput_Mbps_ecdf.mean_x,...
        throughput_Mbps_ecdf.p05);
    % End Faris
end

% --- Executes on button press in constant_axes_limits_checkbox.
function constant_axes_limits_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to constant_axes_limits_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of constant_axes_limits_checkbox
plot_aggregate_results(handles);


% --- Executes on button press in refresh_plots.
function refresh_plots_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_plots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_aggregate_results(handles);


% --- Executes on button press in scatterplot_mean_checkbox.
function scatterplot_mean_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to scatterplot_mean_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of scatterplot_mean_checkbox
plot_aggregate_results(handles);


% --- Executes on button press in scatterplot_points_checkbox.
function scatterplot_points_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to scatterplot_points_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of scatterplot_points_checkbox
plot_aggregate_results(handles);
