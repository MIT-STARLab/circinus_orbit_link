% Author: Kit Kennedy    
%  data rate calculations, as well as outputs for visualization engine

% let's start off by saying this code is terribly structured and hideous. It started out as a small script and I didn't mean to build a bunch of complexity into it, but that has happened at this point. Really should be refactored. I apologize to everyone, now I'm going to go sit in the corner and think about the sadness I've done here :'(

function [rates_output,viz_output] = links_calculator(...
        accesses_data,...
        params,...
        base_directory,...
        verbose)

%  TODO:  at some point should probably add uplink data rate calculation support in here


% TODO: should probably trim off points with data rates equal to zero. For
% now though, just use the safety net provided by the global planner algs,
% where any passes with 0 data vol will be trimmed away

if verbose
    disp('running links calculator');
end

%% Parameters

obs  = accesses_data.obs;
gslink  = accesses_data.gslink;
gsaer  = accesses_data.gsaer;
xlink  = accesses_data.xlink;
xrange  = accesses_data.xrange;

start_time_dt = parse_iso_datestr(params.scenario_start_utc);


xlnk_max_len_secs = params.general_link_params.xlnk_max_len_s;
dlnk_max_len_secs = params.general_link_params.dlnk_max_len_s;
xlnk_storage_decimation = params.general_link_params.xlnk_decimation;
dlnk_storage_decimation = params.general_link_params.dlnk_decimation;

num_sats = params.num_sats;
num_gs = params.num_gs;
num_targets = params.num_targets;

DTR = pi/180;

% built-in comm types
% note that these are NOT USED currently, because the code for calculating them has been removed
comm_types_builtin = containers.Map;  % this is basically a python dict
comm_types_builtin('Optical-general') = 0;
comm_types_builtin('Ka-band-general') = 1;
comm_types_builtin('X-band-general') = 2;
comm_types_builtin('Optical-Sinclair') = 3;
comm_types_builtin('Optical-OCSD') = 4;
comm_types_builtin('UHF-Cadet') = 6;

%% Set up params struct

params_struct.assume_max_datarate = params.general_link_params.assume_max_datarate;  % allows assumption of max data rate on downlink, so it doesn't vary due to rain fade and range.


%% Set up output containers

t_0 = mjuliandate(start_time_dt);

%  output data struct
rates_output = struct;

% Temporary containers
gslink_rates = cell(num_sats,num_gs);
xlink_rates = cell(num_sats,num_sats);
xlink_directions = cell(num_sats,num_sats);

%  data rates outputs
xlink_update = cell(num_sats,num_sats);
xlink_rates_update  = cell(num_sats,num_sats);
xlnkrange_update = cell(num_sats,num_sats);

gslink_update = cell(num_sats,num_gs);
gsaer_update = cell(num_sats,num_gs);
gslink_rates_update = cell(num_sats,num_gs);


%  visualization outputs
t_o = {num_sats,1};
o_locations = {num_sats,1};
t_d = cell(num_sats,1);
d_part = cell(num_sats,1);
t_x = cell(num_sats,1);
x_part = cell(num_sats,1);
dlnk_rate_history = cell(num_sats,num_gs);
xlnk_rate_history = cell(num_sats,num_sats);

%  these structures indicate which downlink and cross-link pairs are to be discluded from processing in the code below. if a particular satellite ID is present as a field within the structure, then the list of disableded ground station or satellite IDs are strings within a cell array for that field.
link_direc_disables_dlnk = params.link_disables.dlnk_direc_disabled_gs_indx_by_sat_indx;
link_direc_disables_xlnk = params.link_disables.xlnk_direc_disabled_xsat_indx_by_sat_indx;

gs_params = params.gs_params;

%% Run

for sat_num = 1:num_sats
    if verbose
        disp(sprintf('sat num: %d',sat_num));
    end

    sat_py_indx_str = ['indx',num2str(sat_num-1)];
    has_dlnk_direc_disables = isfield(link_direc_disables_dlnk,sat_py_indx_str);
    has_xlnk_direc_disables_sat = isfield(link_direc_disables_xlnk,sat_py_indx_str);
    if has_dlnk_direc_disables
        link_direc_disables_dlnk_sat = link_direc_disables_dlnk.(sat_py_indx_str);
    end
    if has_xlnk_direc_disables_sat
        link_direc_disables_xlnk_sat = link_direc_disables_xlnk.(sat_py_indx_str);
    end

    %  get parameters for the satellite

    sat_link_params = params.sat_link_params{sat_num};
    P_tx_dlnk_sat = sat_link_params.dlnk_params.P_tx_W;
    HPBW_dlnk_sat = sat_link_params.dlnk_params.HPBW_rad/1000;
    P_tx_xlnk_sat = sat_link_params.xlnk_params.P_tx_W;
    HPBW_xlnk_sat = sat_link_params.xlnk_params.HPBW_rad/1000;
    pointing_error_sat = sat_link_params.pointing_error_deg * DTR;

    % the comm type for dlnk, specifies which params to use
    comm_type_dlnk = sat_link_params.dlnk_params.comm_type.name;
    % specifies whether or not we're using a built-in model for the downlink link budget or a lookup table from the input params
    comm_type_dlnk_builtin = sat_link_params.dlnk_params.comm_type.built_in;
    comm_type_xlnk_sat = sat_link_params.xlnk_params.comm_type.name;
    comm_type_xlnk_sat_builtin = sat_link_params.xlnk_params.comm_type.built_in;

    [comm_type_dlnk_indx,dlnk_table_range_rates,dlnk_table_ranges,dlnk_table_interpolation] = get_sat_dlnk_params(comm_type_dlnk,comm_type_dlnk_builtin,comm_types_builtin,params);

    [comm_type_xlnk_sat_indx,xlnk_sat_table_range_rates,xlnk_sat_table_ranges,xlnk_sat_table_interpolation] = get_sat_xlnk_params(comm_type_xlnk_sat,comm_type_xlnk_sat_builtin,comm_types_builtin,params);

        
    % Handle observations - note this is only for viz in cesiumJS

    if verbose
        disp('obs');
    end

    obs_count = 1;
    for targ_num = 1:num_targets        
        obs_list = obs{sat_num,targ_num}; 

        num_obs = size(obs_list,1);

        for obs_num = 1:num_obs
            tstart = obs_list(obs_num,1);
            tstop = obs_list(obs_num,2);

            t_o{sat_num,obs_count} =[tstart,tstop];
            o_locations{sat_num,obs_count} = [uint16(targ_num)];

            obs_count = obs_count+1;
        end 
    end

    % have to sort t_o, because obs show times are not created uniquely for each sat-obs target combination, instead they're all lumpted into a single satellite (in VizInputsGenerator.py for cesiumJS input
    % file)
    % Sort by start time
    start_times = zeros(obs_count-1,1);
    for i=1:obs_count-1
        start_times(i) = t_o{sat_num,i}(1);
    end
    [trash , idx] = sort(start_times, 'ascend');

    % now some annoying code to replace t_d for sat_num with sorted
    t_o_new = cell(obs_count-1,1);
    for i=1:obs_count-1
        t_o_new{i} = t_o{sat_num,idx(i)}; % store in new cell array in sorted order
    end
    for i=1:obs_count-1
        t_o{sat_num,i} = t_o_new{i};  % replace back to old cell array
    end
    % end of sorting

    % note: funky things might still happen with display of t_o in cesium
    % if there's any overlap in observation targets...



    %% Handle downlinks

    if verbose
        disp('dlnks');
    end

    % calculate datarates at each timestep during downlinks - store in new cell array gslink_rates
    for gs_num = 1:num_gs
        gs_py_indx_str = num2str(gs_num-1);

        station_params = gs_params.stations{gs_num};
        comm_type_gs =  station_params.comm_type;
        %  check if the dlnk types are compatible ( their names are equal).  if not, skip this gs
        if ~strcmp(comm_type_dlnk,comm_type_gs)
            if verbose
                disp(['skipping gs python indx',gs_py_indx_str,' comm_type does not match'])
            end
            continue
        end

        % link_direc_disables_dlnk_sat
        if has_dlnk_direc_disables
            % check if the ground station index is within the disabled list for the satellite. If yes then skip this ground station
            if any(strcmp(link_direc_disables_dlnk_sat,gs_py_indx_str))
                if verbose
                    disp(['skipping gs python indx',gs_py_indx_str,' comm_type does not match'])
                end
                continue
            end
        end

        dlnks = gslink{sat_num,gs_num};   % gslink is loaded from file above
        dlnks_aer = gsaer{sat_num,gs_num};  % gsaer is loaded from file above

        dlnks_rates = cell(size(dlnks_aer));

        num_dlnks = size(dlnks,1);

        % TODO: add compatibility check between sat and gs comm types

        for dlnk_num = 1:num_dlnks

            dlnk_aer = dlnks_aer{dlnk_num};
            num_timepoints = size(dlnk_aer,1);
            dlnk_rates = zeros(num_timepoints,2);


            for i=1:num_timepoints

                dlnk_rates(i,1) = dlnk_aer(i,1);  % store time
                elev = dlnk_aer(i,3); % in degrees
                range_km = dlnk_aer(i,4);

                % if using built-in link budget
                if comm_type_dlnk_builtin
                    range_m = range_km*1000; % in meters
                    %  note that this immediately gets converted to Mbps
                    dlnk_rates(i,2) = downlink_capacity(P_tx_dlnk_sat,HPBW_dlnk_sat,pointing_error_sat,range_m,comm_type_dlnk_indx,elev,gs_num,params_struct)/1000000;  

                % if using a provided lookup table
                else
                    if strcmp (dlnk_table_interpolation,'floor') 
                        %  get all ranges in the table less than our current range, then take the last of those ("flooring" the current range)
                        rate_entry_indx_unfilt = find(range_km-dlnk_table_ranges>=0);
                        rate_entry_indx = rate_entry_indx_unfilt(end);
                    end

                    dlnk_rates(i,2) = dlnk_table_range_rates(rate_entry_indx,2);
                end

                % leaving in this, emily's old diagnostic code, just in
                % case
                if imag(dlnk_rates(i,2))>0
                    disp('imaginary part in dlnk rate...');
                    range_m
                    elev
                    pause
                end
            end    

            dlnks_rates{dlnk_num} = dlnk_rates;
        end 

        gslink_rates{sat_num,gs_num} = dlnks_rates;
    end

    % TODO: probably ought to trim off any zero downlink data rate time
    % from gslink, gsaer, gslink_rates...left a warning in the code below
    % to catch and warn about this situation...

    % divide downlink windows into smaller time periods to make them more
    % manageable, as well as decimate the gsaer and data rates storage
    dlnk_max_len_dayf = dlnk_max_len_secs/86400;
    for gs_num = 1:num_gs % want to start from sat_num+1 because xlink is symmetric

        gs_py_indx_str = num2str(gs_num-1);

        station_params = gs_params.stations{gs_num};
        comm_type_gs =  station_params.comm_type;
        %  check if the dlnk types are compatible ( their names are equal).  if not, skip this gs
        if ~strcmp(comm_type_dlnk,comm_type_gs)
            continue
        end

        if has_dlnk_direc_disables
            % check if the ground station index is within the disabled list for the satellite. If yes then skip this ground station
            if any(strcmp(link_direc_disables_dlnk_sat,gs_py_indx_str))
                continue
            end
        end

        dlnks = gslink{sat_num,gs_num};   % gslink is loaded from file above
        dlnks_aer = gsaer{sat_num,gs_num};  % gsaer is loaded from file above
        dlnks_rates = gslink_rates{sat_num,gs_num};

        num_dlnks = size(dlnks,1);

        dlnks_update = [];
        dlnks_aer_update = {};
        dlnks_rates_update = {};

        for dlnk_num = 1:num_dlnks

            dlnk = dlnks(dlnk_num,:);
            dlnk_len = dlnk(3);

%             if dlnk_len > dlnk_max_len_secs - % comment out because we
%             want to run this code regardless - it just seemed like the
%             easiest way to include decimation of aer, rates
                % grab these for slicing purposes
                dlnk_aer = dlnks_aer{dlnk_num};
                dlnk_rates = dlnks_rates{dlnk_num};
                dlnk_aer_indx = 1;
                dlnk_aer_indx_start = 1;

                t = dlnk_max_len_dayf; 

                dlnk_start = dlnk(1);
                dlnk_end = dlnk(2);

                last_t = 0;

                % loop through as long as end time of slice is less than
                % the end of the full xlnk access
                while dlnk_start + t < dlnk_end

                    % figure out up to what index in xlnk_xrange (same for
                    % xlnk_rates) we need to slice
                    new_end = dlnk_start + t;
                    dlnk_aer_indx_start = dlnk_aer_indx;
                    while dlnk_aer(dlnk_aer_indx,1) <= new_end                        
                        dlnk_aer_indx = dlnk_aer_indx + 1;
                    end

                    % before we add the new slice, check if there are enough indices for it to be relevant to add
                    if (dlnk_aer_indx-dlnk_aer_indx_start < dlnk_storage_decimation)
                        break
                    end

                    % new sliced access window
                    dlnks_update = [dlnks_update; [dlnk_start + last_t, dlnk_start + t] ];

                    % slice it. Note the -1 so we don't go past the end of
                    % the slice
                    % note the decimation factor used, xlnk_storage_decimation
                    dlnk_aer_slice = dlnk_aer(dlnk_aer_indx_start:dlnk_storage_decimation:dlnk_aer_indx-1,:);
                    dlnk_rates_slice = dlnk_rates(dlnk_aer_indx_start:dlnk_storage_decimation:dlnk_aer_indx-1,:);
                    dlnks_aer_update = [dlnks_aer_update; dlnk_aer_slice];
                    dlnks_rates_update = [dlnks_rates_update; dlnk_rates_slice];

                    last_t = t;
                    t = t + dlnk_max_len_dayf; 
                end

                % handle last slice
                dlnks_update = [dlnks_update; [dlnk_start + last_t, dlnk_end] ];
                dlnk_aer_indx_start = dlnk_aer_indx;
                dlnk_aer_slice = dlnk_aer(dlnk_aer_indx_start:dlnk_storage_decimation:end,:);
                dlnk_rates_slice = dlnk_rates(dlnk_aer_indx_start:dlnk_storage_decimation:end,:);
                dlnks_aer_update = [dlnks_aer_update; dlnk_aer_slice];
                dlnks_rates_update = [dlnks_rates_update; dlnk_rates_slice];

%             else % if we don't need to slice because xlnk isn't too long
%                 dlnks_update = [dlnks_update; dlnk];
%                 dlnks_aer_update = [dlnks_aer_update; dlnks_aer{dlnk_num}];
%                 dlnks_rates_update = [dlnks_rates_update; dlnks_rates{dlnk_num}];
%             end
        end

        gslink_update{sat_num,gs_num} = dlnks_update;
        gsaer_update{sat_num,gs_num} = dlnks_aer_update;
        gslink_rates_update{sat_num,gs_num} = dlnks_rates_update;
    end

    % format dlnk history for input to CesiumJS. Note of course that this
    % is super inefficient, but this code only runs once, and splitting all
    % this code out helps preserve my sanity
    dl_count = 1;
    for gs_num = 1:num_gs

        gs_py_indx_str = num2str(gs_num-1);

        station_params = gs_params.stations{gs_num};
        comm_type_gs =  station_params.comm_type;
        %  check if the dlnk types are compatible ( their names are equal).  if not, skip this gs
        if ~strcmp(comm_type_dlnk,comm_type_gs)
            continue
        end

        if has_dlnk_direc_disables
            % check if the ground station index is within the disabled list for the satellite. If yes then skip this ground station
            if any(strcmp(link_direc_disables_dlnk_sat,gs_py_indx_str))
                continue
            end
        end

        dlnks = gslink_update{sat_num,gs_num};  
        dlnks_aer = gsaer_update{sat_num,gs_num}; 
        dlnks_rates = gslink_rates_update{sat_num,gs_num};

        num_dlnks = size(dlnks,1);

        for dlnk_num = 1:num_dlnks

            dlnk = dlnks(dlnk_num,:);
            dlnk_aer = dlnks_aer{dlnk_num};
            num_timepoints = size(dlnk_aer,1);
            dlnk_rates = dlnks_rates{dlnk_num};

            % store in dlnk history for visualizing in cesiumJS
            tstart = dlnk(1,1); % Need to use dlnk for this point because times in dlnk_aer are decimated
            new_time = (tstart-t_0)*86400 - 0.5;  % 0 rate right before start of downlink 
            dlnk_rate_history{sat_num,gs_num} = [dlnk_rate_history{sat_num,gs_num}; new_time - 0.1,0];   % (minus 0.1 for ~impulsive transition)
            last_history_val = 0;

            for i=1:num_timepoints
                if verbose
                    if dlnk_rates(i,2) == 0
                        % disp('warning, dlnk_rates(i) is zero. Perhaps this point should be trimmed...');
                    end
                end

                % decimate dlnk rate storage history, to maintain some
                % sanity
                % note decimating on top of decimation already done above
                if mod(i-1,2) % note decimating on top of decimation already done above)==0
                    new_time = (dlnk_aer(i,1)-t_0)*86400;
                    dlnk_rate_history{sat_num,gs_num} = [dlnk_rate_history{sat_num,gs_num}; new_time, last_history_val];  % add point for rate we're transitioning from (so transition is impulsive)                            
                    dlnk_rate_history{sat_num,gs_num} = [dlnk_rate_history{sat_num,gs_num}; new_time + 0.1, dlnk_rates(i,2)]; % add point for date rate during pass (plus 0.1 for ~impulsive transition)
                    last_history_val = dlnk_rates(i,2);
                end
            end           

            tstop = dlnk(1,2); % Need to use dlnk for this point because times in dlnk_aer are decimated
            new_time = (tstop-t_0)*86400;
            dlnk_rate_history{sat_num,gs_num} = [dlnk_rate_history{sat_num,gs_num}; new_time + 0.2, last_history_val];  % make sure to add in a point at end of pass, so final rate dropoff is "impulsive". 0.2 so it'll come after any point added above (could possibly conflict with +0.1 above)
            dlnk_rate_history{sat_num,gs_num} = [dlnk_rate_history{sat_num,gs_num}; new_time + 0.3, 0];  % now zero it (plus 0.3 so it'll be, again, after any points above)


            t_d{sat_num,dl_count} =[tstart,tstop];
            % explicitly save partner variables as integers
            d_part{sat_num,dl_count} = [uint16(gs_num)];

            dl_count = dl_count+1;
        end 

    end


    % note: don't need to sort t_d, cesium generation code can
    % handle it, because downlink show times are created separately for
    % each sat-gs combination in VizInputsGenerator.py for cesiumJS input
    % file


    %% Handle crosslinks

    if verbose
        disp('xlnks');
    end

    skip_xlnk = zeros(num_sats,num_sats);
    % calculate datarates at each timestep during crosslinks - store in new cell array xlink_rates
    for xsat_num = sat_num+1:num_sats % want to start from sat_num+1 because xlink is symmetric
        xsat_py_indx_str = num2str(xsat_num-1);

        xsat_link_params = params.sat_link_params{xsat_num};
        comm_type_xlnk_xsat = xsat_link_params.xlnk_params.comm_type.name;

        skip_xlnk_sat_xsat = false;
        skip_xlnk_xsat_sat = false;

        %  check if the crosslink types are compatible ( their names are equal).  if not, skip this satellite
        % TODO: If different technologies are desired for TX and RX, need to update this check here. for now are assuming that both satellites must have the same technology
        if ~strcmp(comm_type_xlnk_sat,comm_type_xlnk_xsat)
            if verbose
                disp(['skipping xsat python indx',xsat_py_indx_str])
            end
            continue
        end

        %  determine if the cross-link is explicitly disabled in either direction. if yes, then we'll be skipping a direction
        has_xlnk_direc_disables_xsat = isfield(link_direc_disables_xlnk,xsat_py_indx_str);
        if has_xlnk_direc_disables_xsat
            link_direc_disables_xlnk_xsat = link_direc_disables_xlnk.(xsat_py_indx_str);
        end

        if has_xlnk_direc_disables_sat
            % check if the xsat index is within the disabled list for the satellite. If yes then skip this satellite
            skip_xlnk_sat_xsat = skip_xlnk_sat_xsat || any(strcmp(link_direc_disables_xlnk_sat,xsat_py_indx_str));
            
        end

        if has_xlnk_direc_disables_xsat
            % check if the xsat index is within the disabled list for the satellite. If yes then skip this satellite
            skip_xlnk_xsat_sat = skip_xlnk_xsat_sat || any(strcmp(link_direc_disables_xlnk_xsat,sat_py_indx_str));
            
        end

        if skip_xlnk_sat_xsat
            if verbose
                disp(['skipping xlnk ',sat_py_indx_str,'->',xsat_py_indx_str,' (py indcs)'])
            end
        end
        if skip_xlnk_xsat_sat
            if verbose
                disp(['skipping xlnk ',xsat_py_indx_str,'->',sat_py_indx_str,' (py indcs)'])
            end
        end

        %  only don't procede through the rest of this iteration  if we are skipping both directions in the cross-link
        if skip_xlnk_sat_xsat && skip_xlnk_xsat_sat
            skip_xlnk(sat_num,xsat_num) = 1;
            continue
        end


        P_tx_xlnk_xsat = xsat_link_params.xlnk_params.P_tx_W;
        HPBW_xlnk_xsat = xsat_link_params.xlnk_params.HPBW_rad/1000;
        pointing_error_xsat = xsat_link_params.pointing_error_deg * DTR;

        comm_type_xlnk_xsat_builtin = xsat_link_params.xlnk_params.comm_type.built_in;
        [comm_type_xlnk_xsat_indx,xlnk_xsat_table_range_rates,xlnk_xsat_table_ranges,xlnk_xsat_table_interpolation] = get_sat_xlnk_params(comm_type_xlnk_xsat,comm_type_xlnk_xsat_builtin,comm_types_builtin,params);

        %  record if they are the same params, in which case the calculated data rates will be the same so we can save some processing
        same_params_both_dirs = (P_tx_xlnk_sat==P_tx_xlnk_xsat) && (HPBW_xlnk_sat==HPBW_xlnk_xsat) && (pointing_error_sat==pointing_error_xsat) && strcmp(comm_type_xlnk_sat,comm_type_xlnk_xsat);

        xlnks = xlink{sat_num,xsat_num};   % xlink is loaded from file above
        xlnks_xrange = xrange{sat_num,xsat_num};  % xrange is loaded from file above

        xlnks_rates = cell(size(xlnks_xrange));
        xlnks_direcs = cell(size(xlnks_xrange));

        num_xlnks = size(xlnks,1);

        for xlnk_num = 1:num_xlnks

            xlnk_xrange = xlnks_xrange{xlnk_num};
            num_timepoints = size(xlnk_xrange,1);
            %  for each row, first index is timestamp, second index is data rate from satellite to cross-link satellite, and third index is data rate from cross-link satellite to satellite
            xlnk_rates = zeros(num_timepoints,3);

            for i=1:num_timepoints
                xlnk_rates(i,1) = xlnk_xrange(i,1);

                range_km = xlnk_xrange(i,2);

                if ~comm_type_xlnk_sat_builtin && strcmp (xlnk_sat_table_interpolation,'floor') 
                    %  get all ranges in the table less than our current range, then take the last of those ("flooring" the current range)
                    rate_entry_indx_unfilt = find(range_km-xlnk_sat_table_ranges>=0);
                    rate_entry_indx_sat = rate_entry_indx_unfilt(end);
                end
                
                if ~comm_type_xlnk_xsat_builtin && strcmp (xlnk_xsat_table_interpolation,'floor') 
                    %  get all ranges in the table less than our current range, then take the last of those ("flooring" the current range)
                    rate_entry_indx_unfilt = find(range_km-xlnk_xsat_table_ranges>=0);
                    rate_entry_indx_xsat = rate_entry_indx_unfilt(end);
                end

                if skip_xlnk_sat_xsat
                    xlnk_rates(i,2) = 0;
                else
                    if comm_type_xlnk_sat_builtin
                        xlnk_rates(i,2) = calc_xlnk_rate(range_km,P_tx_xlnk_sat,HPBW_xlnk_sat,pointing_error_sat,range_m,comm_type_xlnk_sat_indx,params_struct,comm_type_xlnk_sat_builtin,xlnk_sat_table_interpolation,xlnk_sat_table_ranges,xlnk_sat_table_range_rates);
                    else
                        xlnk_rates(i,2) = xlnk_sat_table_range_rates(rate_entry_indx_sat,2);
                    end
                end

                % FML, this is so ugly...
                if skip_xlnk_xsat_sat
                    xlnk_rates(i,3) = 0;
                else
                    if same_params_both_dirs && ~skip_xlnk_sat_xsat
                        xlnk_rates(i,3) = xlnk_rates(i,2);
                    else
                        if comm_type_xlnk_xsat_builtin
                            xlnk_rates(i,3) = calc_xlnk_rate(range_km,P_tx_xlnk_xsat,HPBW_xlnk_xsat,pointing_error_xsat,range_m,comm_type_xlnk_xsat_indx,params_struct,comm_type_xlnk_xsat_builtin,xlnk_xsat_table_interpolation,xlnk_xsat_table_ranges,xlnk_xsat_table_range_rates);
                        else
                            xlnk_rates(i,3) = xlnk_xsat_table_range_rates(rate_entry_indx_xsat,2);
                        end
                    end
                end

                % leaving in this, emily's old diagnostic code, just in
                % case
                if imag(xlnk_rates(i,2))>0
                    disp('imaginary part in xlnk rate...setting to 73 mbps for now');

                    xlnk_rates(i,2) = 73.1133803921651;
                end
                if imag(xlnk_rates(i,3))>0
                    disp('imaginary part in xlnk rate...setting to 73 mbps for now');

                    xlnk_rates(i,3) = 73.1133803921651;
                end
            end    

            xlnk_direcs = zeros(1,3);
            if ~skip_xlnk_sat_xsat
                %  record a "true" for the first satellite transmitting
                xlnk_direcs(1,1) = 1;
            end
            if ~skip_xlnk_xsat_sat
                %  record a "true" for the second satellite transmitting
                xlnk_direcs(1,2) = 1;
            end
            if same_params_both_dirs
                %  record a "true" for the satellites transmitting at the same rate
                xlnk_direcs(1,3) = 1;
            end

            xlnks_rates{xlnk_num} = xlnk_rates;
            xlnks_direcs{xlnk_num} = xlnk_direcs;
        end 

        xlink_rates{sat_num,xsat_num} = xlnks_rates;
        xlink_directions{sat_num,xsat_num} = xlnks_direcs;
    end


    % divide crosslink windows into smaller time periods to make them more
    % manageable, as well as decimate the xrange and data rates storage
    xlnk_max_len_dayf = xlnk_max_len_secs/86400;
    for xsat_num = sat_num+1:num_sats % want to start from sat_num+1 because xlink is symmetric
        
        if skip_xlnk(sat_num,xsat_num) == true
            continue
        end
        
        xlnks = xlink{sat_num,xsat_num};   % xlink is loaded from file above
        xlnks_xrange = xrange{sat_num,xsat_num};  % xrange is loaded from file above
        xlnks_rates = xlink_rates{sat_num,xsat_num};
        xlnks_direcs = xlink_directions{sat_num,xsat_num};

        num_xlnks = size(xlnks,1);

        xlnks_update = [];
        xlnks_xrange_update = {};
        xlnks_rates_update = {};

        for xlnk_num = 1:num_xlnks

            xlnk = xlnks(xlnk_num,:);
            xlnk_len = xlnk(3);

%             if xlnk_len > xlnk_max_len_secs  % comment out because we
%             want to run this code regardless - it just seemed like the
%             easiest way to include decimation of xrange, rates
                % grab these for slicing purposes
                xlnk_xrange = xlnks_xrange{xlnk_num};
                xlnk_rates = xlnks_rates{xlnk_num};
                xlnk_direcs = xlnks_direcs{xlnk_num};
                xlnk_xrange_indx = 1;
                xlnk_xrange_indx_start = 1;

                t = xlnk_max_len_dayf; 

                xlnk_start = xlnk(1);
                xlnk_end = xlnk(2);

                last_t = 0;

                % loop through as long as end time of slice is less than
                % the end of the full xlnk access
                while xlnk_start + t < xlnk_end

                    % figure out up to what index in xlnk_xrange (same for
                    % xlnk_rates) we need to slice
                    new_end = xlnk_start + t;
                    xlnk_xrange_indx_start = xlnk_xrange_indx;
                    while xlnk_xrange(xlnk_xrange_indx,1) <= new_end                        
                        xlnk_xrange_indx = xlnk_xrange_indx + 1;
                    end

                    % before we add the new slice, check if there are enough indices for it to be relevant to add
                    if (xlnk_xrange_indx-xlnk_xrange_indx_start < xlnk_storage_decimation)
                        break
                    end

                    % new sliced access window
                    % combine with the cross-linked directions array
                    xlnks_update = [xlnks_update; [[xlnk_start + last_t, xlnk_start + t], xlnk_direcs] ];

                    % slice it. Note the -1 so we don't go past the end of
                    % the slice
                    % note the decimation factor used, xlnk_storage_decimation
                    xlnk_xrange_slice = xlnk_xrange(xlnk_xrange_indx_start:xlnk_storage_decimation:xlnk_xrange_indx-1,:);
                    xlnk_rates_slice = xlnk_rates(xlnk_xrange_indx_start:xlnk_storage_decimation:xlnk_xrange_indx-1,:);
                    
                    xlnks_xrange_update = [xlnks_xrange_update; xlnk_xrange_slice];
                    xlnks_rates_update = [xlnks_rates_update; xlnk_rates_slice];

                    last_t = t;
                    t = t + xlnk_max_len_dayf; 
                end

                % handle last slice
                xlnks_update = [xlnks_update; [[xlnk_start + last_t, xlnk_end], xlnk_direcs] ];
                xlnk_xrange_indx_start = xlnk_xrange_indx;
                xlnk_xrange_slice = xlnk_xrange(xlnk_xrange_indx_start:xlnk_storage_decimation:end,:);
                xlnk_rates_slice = xlnk_rates(xlnk_xrange_indx_start:xlnk_storage_decimation:end,:);
                xlnks_xrange_update = [xlnks_xrange_update; xlnk_xrange_slice];
                xlnks_rates_update = [xlnks_rates_update; xlnk_rates_slice];

%             else % if we don't need to slice because xlnk isn't too long
%                 xlnks_update = [xlnks_update; xlnk];
%                 xlnks_xrange_update = [xlnks_xrange_update; xlnks_xrange{xlnk_num}];
%                 xlnks_rates_update = [xlnks_rates_update; xlnks_rates{xlnk_num}];
%             end
        end

        xlink_update{sat_num,xsat_num} = xlnks_update;
        xlnkrange_update{sat_num,xsat_num} = xlnks_xrange_update;
        xlink_rates_update{sat_num,xsat_num} = xlnks_rates_update;
    end

    % TODO: this code currently ignores the fact the update to cross-link rates that  allows the two directions between satellitesto be nonsymmetric. should update this to include it

    % format xlnk history for input to CesiumJS. Note of course that this
    % is super inefficient, but this code only runs once, and splitting all
    % this code out helps preserve my sanity
    xl_count = 1;
    for xsat_num = sat_num+1:num_sats % want to start from sat_num+1 because xlink is symmetric
        if skip_xlnk(sat_num,xsat_num) == true
            continue
        end

        xlnks = xlink_update{sat_num,xsat_num};   
        xlnks_xrange = xlnkrange_update{sat_num,xsat_num}; 
        xlnks_rates = xlink_rates_update{sat_num,xsat_num};

        num_xlnks = size(xlnks,1);

        for xlnk_num = 1:num_xlnks

            xlnk= xlnks(xlnk_num,:);
            xlnk_xrange = xlnks_xrange{xlnk_num};
            num_timepoints = size(xlnk_xrange,1);
            xlnk_rates = xlnks_rates{xlnk_num};

            % store in xlnk history for visualizing in cesiumJS
            tstart = xlnk(1,1); % Need to use xlnk for this point because times in xlnk_xrange are decimated
            new_time = (tstart-t_0)*86400 - 0.5;  % 0 rate right before start of crosslink 
            xlnk_rate_history{sat_num,xsat_num} = [xlnk_rate_history{sat_num,xsat_num}; new_time - 0.1,0];   % (minus 0.1 for ~impulsive transition)
            last_history_val = 0;

            for i=1:num_timepoints
                if verbose
                    if xlnk_rates(i,2) == 0
                        % disp('warning, xlnk_rates(i,2) is zero. Perhaps this point should be trimmed...');
                    end
                end

                % decimate xlnk rate storage history, to maintain some
                % sanity
                % note decimating on top of decimation already done above
                if mod(i-1,2)==0
                    new_time = (xlnk_xrange(i,1)-t_0)*86400; % use time from range data struct because that captures every point in the window
                    xlnk_rate_history{sat_num,xsat_num} = [xlnk_rate_history{sat_num,xsat_num}; new_time, last_history_val];  % add point for rate we're transitioning from (so transition is impulsive)                            
                    xlnk_rate_history{sat_num,xsat_num} = [xlnk_rate_history{sat_num,xsat_num}; new_time + 0.1, xlnk_rates(i,2)]; % add point for date rate during pass (plus 0.1 for ~impulsive transition)
                    last_history_val = xlnk_rates(i,2);
                end
            end           

            tstop = xlnk(1,2); % Need to use xlnk for this point because times in xlnk_xrange are decimated
            new_time = (tstop-t_0)*86400;
            xlnk_rate_history{sat_num,xsat_num} = [xlnk_rate_history{sat_num,xsat_num}; new_time + 0.2, last_history_val];  % make sure to add in a point at end of pass, so final rate dropoff is "impulsive". 0.2 so it'll come after any point added above (could possibly conflict with +0.1 above)
            xlnk_rate_history{sat_num,xsat_num} = [xlnk_rate_history{sat_num,xsat_num}; new_time + 0.3, 0];  % now zero it (plus 0.3 so it'll be, again, after any points above)


            t_x{sat_num,xl_count} =[tstart,tstop];
            % explicitly save partner variables as integers
            x_part{sat_num,xl_count} = [uint16(xsat_num)];

            xl_count = xl_count+1;
        end 

    end

    % note: don't need to sort t_x, cesium generation code can
    % handle it, because downlink show times are created separately for
    % each sat-gs combination in VizInputsGenerator.py for cesiumJS input
    % file

end

%% Save output
rates_output.obs = obs;
rates_output.xlink_update = xlink_update;
rates_output.xlink_rates_update = xlink_rates_update;
rates_output.xlnkrange_update = xlnkrange_update;
rates_output.gslink_update = gslink_update;
rates_output.gsaer_update = gsaer_update;
rates_output.gslink_rates_update = gslink_rates_update;


% save outputs for cesiumJS viz
viz_output.t_o =t_o;
viz_output.o_locations =o_locations;
viz_output.t_d =t_d;
viz_output.d_part =d_part;
viz_output.t_x =t_x;
viz_output.x_part =x_part;
viz_output.dlnk_rate_history =dlnk_rate_history;
viz_output.xlnk_rate_history =xlnk_rate_history;

% q_o_sizes_history
% gs_availability_windows
% batt_stored_history
% t_eclipse

