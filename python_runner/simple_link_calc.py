
import numpy as np
from matlab_links_wrapper import matlab_links_wrapper

# data_rates file setup, specific numbers correspond to zhou use case
""" data_rates_example['accesses_data_rates'] - dict of dicts, relevant ones are:
'dlnk_rates':
- list of N_sat=6 lists, where each list is:
-- [0] list of N_ground_stations=4 lists, where each list is:
--- [0][0] list of N_access_windows=4 lists, where each list is:
---- [0][0][0] list of N_timesteps_per_window lists, where each list is:
----- [0][0][0][0] list with two elements (MJD_timestep, data_rate_Mbps)
'dlnk_times':
- list of N_sat=6 lists, where each list is:
-- [0] list of N_ground_stations=4 lists, where each list is:
--- [0][0] list of N_access_windows=4 lists, where each list is:
---- [0][0][0] list with two elements (JD_start,JD_end)
'obs_times': (same as orbit_prop_data['accesses_data']['obs_times'])
- list of N_sat=6 lists, where each list is:
-- [0] list of N_obs_targets=4 lists, where each list is:
--- [0][3] list of N_access_windows=2 lists, where each list is:
---- [0][3][0] list with two elements (JD_start,JD_end, duration_seconds)
'xlnk_rates':
- list of N_sat=6 lists, where each list is:
-- [0] list of N_xsat=6 lists, where each list is:
--- [0][3] list of N_access_windows=432 lists, where each list is:
---- [0][3][0] list of N_timesteps_per_window=4 lists, where each list is:
----- [0][3][0][0] list with three elements (MJD_timestep, data_rate_Mbps sat->xsat, data_rate_Mbps xsat->sat)
'xlnk_times':
- list of N_sat=6 lists, where each list is:
-- [0] list of N_xsat=6 lists, where each list is:
--- [0][3] list of N_access_windows=432 lists, where each list is:
---- [0][3][0] list with five elements (MJD_start,MJD_end, bool_1, bool_2, bool_3), where
----- bool_1: true if the first satellite (sat) is transmitting
----- bool_2: true if the second satellite (xsat) is transmitting
----- bool_3: in the case that they're both transmitting, true if they're transmitting at the same rate over the entire window. Note that this can also be true if only one is transmitting, so it's a necessary but not sufficient condition for symmetry"
"""

SEC_TO_DAY = 1/(24*60*60)

def pack_dlnks(orbit_prop_dlnk_times,dlnk_rate_Mbps):
    '''pack up downlinks into correct format based on max dlnk rate rule'''
    
    window_step_seconds = 20
    # these time between each MJD in the "timesteps" is 20 seconds, not 10 seconds, and it correspponds to the end of the pair of timesteps at the lower level

    dlnk_times = []
    dlnk_rates = []
    for sat_ind in range(len(orbit_prop_dlnk_times)):
        dlnk_rates.append([])
        dlnk_times.append([])
        for gs_ind in range(len(orbit_prop_dlnk_times[sat_ind])):
            dlnk_rates[sat_ind].append([])
            dlnk_times[sat_ind].append([])
            for window_ind in range(len(orbit_prop_dlnk_times[sat_ind][gs_ind])):
                dlnk_rates[sat_ind][gs_ind].append([])
                dlnk_times[sat_ind][gs_ind].append([])
                # the ONLY difference between "dlnk_times" output of orbit_link_public and the orbit_prop_data['accesses_data']['dlnk_times'] is that the duration_seconds (last element) is removed
                dlnk_times[sat_ind][gs_ind][window_ind] = orbit_prop_dlnk_times[sat_ind][gs_ind][window_ind][:-1]
                duration_s = orbit_prop_dlnk_times[sat_ind][gs_ind][window_ind][-1]
                num_window_steps = int(np.floor(duration_s/window_step_seconds))
                MJD_start = orbit_prop_dlnk_times[sat_ind][gs_ind][window_ind][0]
                for step_ind in range(num_window_steps):
                    # the cur_MJD saved seems to the "at the end of each window step"
                    cur_MJD = MJD_start + window_step_seconds*(step_ind+1)*SEC_TO_DAY 
                    dlnk_rates[sat_ind][gs_ind][window_ind].append([cur_MJD, dlnk_rate_Mbps])

    
    return (dlnk_times, dlnk_rates) 
    
def pack_xlnks(orbit_prop_xlnk_times,xlnk_rate_Mbps,xlnk_max_len_s):
    '''pack up downlinks into correct format based on max xlnk rate rule'''

    bool_list = [1, 1, 1]
    """
    ----- bool_1: true if the first satellite (sat) is transmitting
    ----- bool_2: true if the second satellite (xsat) is transmitting
    ----- bool_3: in the case that they're both transmitting, true if they're transmitting at the same rate over the entire window. Note that this can also be true if only one is transmitting, so it's a necessary but not sufficient condition for symmetry"
    """
    window_step_seconds = 50
    # these time between each MJD in the "timesteps" is 50 seconds, not 10 seconds, and it correspponds to the end of the pair of timesteps at the lower level

    xlnk_times = []
    xlnk_rates = []
    for sat_ind in range(len(orbit_prop_xlnk_times)):
        xlnk_rates.append([])
        xlnk_times.append([])
        for xsat_ind in range(len(orbit_prop_xlnk_times[sat_ind])):
            xlnk_rates[sat_ind].append([])
            xlnk_times[sat_ind].append([])
            cur_new_window_ind = 0
            for window_ind in range(len(orbit_prop_xlnk_times[sat_ind][xsat_ind])):
                MJD_start = orbit_prop_xlnk_times[sat_ind][xsat_ind][window_ind][0]
                duration_s = orbit_prop_xlnk_times[sat_ind][xsat_ind][window_ind][2]
                if duration_s <= xlnk_max_len_s:
                    # still only one window
                    xlnk_rates[sat_ind][xsat_ind].append([])
                    xlnk_times[sat_ind][xsat_ind].append([])
                    xlnk_times[sat_ind][xsat_ind][cur_new_window_ind] = [MJD_start, MJD_end] + bool_list
                    window_steps = int(np.floor(duration_s/window_step_seconds))
                    for step_ind in range(window_steps):
                        # the cur_MJD saved seems to the "at the end of each window step"
                        cur_MJD = MJD_start + window_step_seconds*(step_ind+1)*SEC_TO_DAY 
                        xlnk_rates[sat_ind][xsat_ind][cur_new_window_ind].append([cur_MJD,xlnk_rate_Mbps ,xlnk_rate_Mbps])
                    cur_new_window_ind += 1
                else:
                    # need to break the window up into xlnk_max_len_s size windows
                    num_new_windows = int(np.floor(duration_s/xlnk_max_len_s))
                    # NOTE: this method may drop a fraction of a window (< xlnk_max_len_s) at the end
                    # need a to count the previous ending new_window ind and add it:

                    for new_window_ind in range(num_new_windows):
                        xlnk_rates[sat_ind][xsat_ind].append([])
                        xlnk_times[sat_ind][xsat_ind].append([])
                        MJD_start += xlnk_max_len_s*(new_window_ind)*SEC_TO_DAY 
                        MJD_end = MJD_start + xlnk_max_len_s*SEC_TO_DAY 
                        xlnk_times[sat_ind][xsat_ind][cur_new_window_ind] = [MJD_start, MJD_end] + bool_list
                        new_window_steps = int(np.floor(xlnk_max_len_s/window_step_seconds))
                        for step_ind in range(new_window_steps):
                            # the cur_MJD saved seems to the "at the end of each window step"
                            cur_MJD = MJD_start + window_step_seconds*(step_ind+1)*SEC_TO_DAY 
                            xlnk_rates[sat_ind][xsat_ind][cur_new_window_ind].append([cur_MJD,xlnk_rate_Mbps ,xlnk_rate_Mbps])
                        cur_new_window_ind += 1 # for indexing into total array
    return (xlnk_times, xlnk_rates)               


def pack_links_data(orbit_prop_data,dlnk_rate_Mbps,xlnk_rate_Mbps,xlnk_max_len_s):
    '''wrapper to be able to test just data rate packing outside of pipeline'''

    accesses_data_rates = {}

    (dlnk_times, dlnk_rates) = pack_dlnks(orbit_prop_data['accesses_data']['dlnk_times'],dlnk_rate_Mbps)
    (xlnk_times, xlnk_rates) = pack_xlnks(orbit_prop_data['accesses_data']['xlnk_times'],xlnk_rate_Mbps,xlnk_max_len_s)

    # POPULATE ACCESSES DATA
    accesses_data_rates["obs_times"]  = orbit_prop_data['accesses_data']['obs_times']
    accesses_data_rates["xlnk_times"]  = xlnk_times
    accesses_data_rates["xlnk_rates"] = xlnk_rates

    accesses_data_rates["dlnk_times"] = dlnk_times
    accesses_data_rates["dlnk_rates"] = dlnk_rates

    return accesses_data_rates


def py_process_links(sim_case_config, constellation_config,gs_network_config,ops_profile_config,sim_general_config,sat_config,orbit_prop_data,verbose=True):
    '''another wrapper to format outputs'''

    include_xlnk_range_in_output = sim_general_config["general_sim_params"]["general_link_params"]['include_xlnk_range_in_output']
    include_dlnk_aer_in_output = sim_general_config["general_sim_params"]["general_link_params"]['include_dlnk_aer_in_output']

    
    viz_data = {}
    # TODO: pull this from the right spot (not sure where, one of the params)
    dlnk_rate_Mbps = 20
    xlnk_rate_Mbps = 10
    
    # MAX duration for crosslinks and downlinks (TODO: enforce for downlinks, only enforce Xlnks now)
    xlnk_max_len_s = sim_general_config['general_sim_params']['general_link_params']['xlnk_max_len_s']
    #dlnk_max_len_s = sim_general_config['general_sim_params']['general_link_params']['dlnk_max_len_s']

    assume_max_datarate_flag = sim_general_config['general_sim_params']['general_link_params']['assume_max_datarate']

    assume_max_datarate_flag = False

    if assume_max_datarate_flag:
        accesses_data_rates = pack_links_data(orbit_prop_data,dlnk_rate_Mbps,xlnk_rate_Mbps,xlnk_max_len_s)


        accesses_data_rates["_comment_xlnk_times"]  = "First nesting level is each sat index. Second nesting level is each other sat that can be a cross-link partner. Third nesting level is  each individual cross-link within that sat,xsat pair. For each of the 5-element arrays at that third nesting level,  the first double is the start time of the window (in modified Julian date), the second double is the end time of the window, and the remaining three doubles are intended to be logical values: 1. true if the first satellite is transmitting, 2. true if the second satellite is transmitting, 3. and in the case that they're both transmitting, true if they're transmitting at the same rate over the entire window. Note that this can also be true if only one is transmitting, so it's a necessary but not sufficient condition for symmetry"
        accesses_data_rates["_comment_xlnk_rates"]  = "First nesting level is each sat index. Second nesting level is each other sat that can be a cross-link partner (xsat index). Third nesting level is each individual cross-link within that sat,xsat pair. Fourth nesting level is each time point within the cross-link. For each time point there are three fields: 1.  the time (modified Julian date), 2. the transmit data rate in direction sat index -> xsat index, 3. the transmit data rate in direction xsat index -> sat index"
    else:
        accesses_data_rates, viz_data = matlab_links_wrapper(sim_case_config, constellation_config,gs_network_config,ops_profile_config,sim_general_config,sat_config,orbit_prop_data,verbose)
        return accesses_data_rates, viz_data


    if include_xlnk_range_in_output:
        accesses_data_rates["xlnk_range"] = orbit_prop_data['accesses_data']['xlnk_range']
    if include_dlnk_aer_in_output:
        accesses_data_rates["dlnk_aer"] = orbit_prop_data['accesses_data']['dlnk_aer']
    
        
    # POPULATE VIZ DATA (EMPTY FOR NOT SINCE WE ARE NOT USING VISUALIZATION)
    viz_data["obs_times_flat"]  = []
    viz_data["dlnk_times_flat"] = []
    viz_data["xlnk_times_flat"] = []
    viz_data["dlnk_rate_history"] = []
    viz_data["dlnk_rate_history_epoch"] = sim_case_config['scenario_params']['start_utc'] 
    viz_data["xlnk_rate_history"] = []
    viz_data["xlnk_rate_history_epoch"] = sim_case_config['scenario_params']['start_utc'] 

    viz_data["obs_locations"] = []
    viz_data["dlnk_partners"] = []
    viz_data["xlnk_partners"] = []

    return accesses_data_rates, viz_data

def py_links_wrapper(data):
    ''' links wrapper to match interface and outputs form'''
    # TODO: there is some overlap here with my code

    orbit_prop_data = data['orbit_prop_data']
    sim_case_config = data['sim_case_config'] # TODO; why pack it into data and then re-unpack it once in run? 
    constellation_config = data['constellation_config']
    gs_network_config = data['gs_network_config']
    ops_profile_config = data['ops_profile_config']
    sim_general_config = data['sim_general_config']
    sat_config = data['sat_config']

    data_rates_output = {}
    data_rates_output['version'] = "0.3"
    data_rates_output['scenario_params'] = sim_case_config['scenario_params'] 

    viz_output = {}
    viz_output['version'] = "0.1"
    viz_output['scenario_params'] = sim_case_config['scenario_params'] 

    out_stuff = py_process_links(
        sim_case_config,
        constellation_config,
        gs_network_config,
        ops_profile_config,
        sim_general_config,
        sat_config,
        orbit_prop_data,
        verbose=True)

    data_rates_output['accesses_data_rates'] = out_stuff[0]
    # 'other data' only has eclipse_times in it
    data_rates_output['other_data']= {'eclipse_times': orbit_prop_data['accesses_data']['ecl_times']}
    viz_output['viz_data'] = out_stuff[1]

    return data_rates_output, viz_output



# for sat_link_history, may not be necessary and only for visualization purposes.
