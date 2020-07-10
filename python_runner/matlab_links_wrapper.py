""" Does the pseudo-matlab-ifying from process_links in runner_orbitlink.py """
from runner_orbitlink import DATA_RATE_OUTPUT_JSON_VER, VIZ_OUTPUT_JSON_VER
import io_tools
import sys
import links_wrapper


def wrap_dict_arr(d):
    out = {}
    for k, v in d.items():
        out[k] = [v]
    return out


def configure_inputs(sim_case_config, constellation_config,gs_network_config,ops_profile_config,sim_general_config,sat_config,orbit_prop_data,verbose):
    if verbose:
        print('Converting inputs to matlab format')

    num_sats = constellation_config['constellation_definition']['constellation_params']['num_satellites']

    # Prepare sat_link_params
    sat_id_prefix= constellation_config['constellation_definition']['constellation_params']['sat_id_prefix']

    slp_holder = {
        "sat_id_prefix": constellation_config['constellation_definition']['constellation_params']['sat_id_prefix'],
        "sat_ids": constellation_config['constellation_definition']['constellation_params']['sat_ids'],
        **sat_config['sat_model_definition']['sat_model_params']['sat_link_params'][0]
    }

    sat_config['sat_model_definition']['sat_model_params']['sat_link_params'] = [slp_holder]

    sat_id_order = constellation_config['constellation_definition']['constellation_params']['sat_id_order']

    sat_link_params, all_sat_ids = io_tools.unpack_sat_entry_list(
        sat_config['sat_model_definition']['sat_model_params']["sat_link_params"],
        force_duplicate = True)
    sat_id_order = io_tools.make_and_validate_sat_id_order(sat_id_order,sat_id_prefix,num_sats,all_sat_ids)

    io_tools.validate_ids(validator=sat_id_order,validatee=all_sat_ids)

    sat_link_params = io_tools.sort_input_params_by_sat_IDs(sat_link_params,sat_id_order)

    # Helper functions
    # (This has to come after sat link param definition because of sat_id_order)
    def get_sat_indx_str(sat_ID):
        return str(sat_id_order.index(str(sat_ID)))

    def get_gs_indx_str(gs_ID):
        gs_id_order= io_tools.make_and_validate_gs_id_order(gs_network_config['network_definition']['gs_net_params'])
    
        return str(gs_id_order.index(str(gs_ID)))
    
    def convert_list_elems(lst,f_conv):
        return [f_conv(elem) for elem in lst]

    # Prepare sub structures
    # These stringified dictionaries can almost certainly be simplified since it
    # 's not going to matlab
    dlnk_direc_disabled_gs_ID_by_sat_IDstr = ops_profile_config['ops_profile_params']["link_disables"]['dlnk_direc_disabled_gs_ID_by_sat_IDstr']
    xlnk_direc_disabled_xsat_ID_by_sat_IDstr = ops_profile_config['ops_profile_params']["link_disables"]['xlnk_direc_disabled_xsat_ID_by_sat_IDstr']

    dlnk_direc_disabled_gs_indx_by_sat_indx = {'indx'+get_sat_indx_str(sat_ID): convert_list_elems(disable_IDs,get_gs_indx_str) for sat_ID,disable_IDs in dlnk_direc_disabled_gs_ID_by_sat_IDstr.items()}
    
    xlnk_direc_disabled_xsat_indx_by_sat_indx = {'indx'+get_sat_indx_str(sat_ID): convert_list_elems(disable_IDs,get_sat_indx_str) for sat_ID,disable_IDs in xlnk_direc_disabled_xsat_ID_by_sat_IDstr.items()}

    link_disables = {
        "dlnk_direc_disabled_gs_indx_by_sat_indx": dlnk_direc_disabled_gs_indx_by_sat_indx,
        "xlnk_direc_disabled_xsat_indx_by_sat_indx": xlnk_direc_disabled_xsat_indx_by_sat_indx,
    }

    # Actual inputs
    # TODO: reshape according to matlab_nesting settings
    accesses_data = {
        'obs': orbit_prop_data["accesses_data"]["obs_times"],
        'gslink': orbit_prop_data["accesses_data"]["dlnk_times"],
        'gsaer': orbit_prop_data["accesses_data"]["dlnk_aer"],
        'xlink': orbit_prop_data["accesses_data"]["xlnk_times"],
        'xrange': orbit_prop_data["accesses_data"]["xlnk_range"],
    }

    params = { # Coppied from matlab-ify the args section
        'scenario_start_utc': sim_case_config['scenario_params']['start_utc'],
        'num_sats': num_sats,
        'num_gs': 
            gs_network_config['network_definition']['gs_net_params']['num_stations'],
        'num_targets': 
            ops_profile_config['ops_profile_params']['obs_params']['num_targets'],
        'gs_params': gs_network_config['network_definition']['gs_net_params'],
        'general_link_params': 
            wrap_dict_arr(sim_general_config["general_sim_params"]["general_link_params"]),
        'lookup_params': sim_case_config['scenario_params']["lookup_params"],
        'assume_max_datarate': False, # If this is true, we just use simple link
        'link_disables': link_disables,
        'sat_link_params': sat_link_params,
        'verbose': sim_general_config["general_sim_params"][
            "general_link_params"]['matlab_verbose_links'],
    }

    return accesses_data, params


def check_inputs(sim_case_config, constellation_config,gs_network_config,ops_profile_config,sim_general_config,sat_config,orbit_prop_data,verbose):
    if sat_config['sat_model_definition']['sat_model_params']["sat_link_params"][0]["dlnk_params"]["comm_type"]["built_in"]:
        raise NotImplementedError('Currently the code for built in link models is not included in this repository')
    if sat_config['sat_model_definition']['sat_model_params']["sat_link_params"][0]["xlnk_params"]["comm_type"]["built_in"]:
        raise NotImplementedError('Currently the code for built in link models is not included in this repository')
    
    # Not doing version checks the same way anymore.
    # if orbit_prop_data['version'] != "0.3":
    #     # No other versions had 
    #     raise NotImplementedError


def configure_outputs(rates_output_by_sat_ml, viz_output_by_sat_ml, sim_case_config, constellation_config,gs_network_config,ops_profile_config,sim_general_config,sat_config,orbit_prop_data,verbose):
    """
    rates_output_by_sat_ml, viz_output_by_sat_ml: outputs from output matlab func

    sim_case_config, constellation_config,gs_network_config,ops_profile_config,sim_general_config,sat_config,orbit_prop_data,verbose: inputs to original function
    """

    # This part is taken from PipelineRunner run function
    data_rates_output = {}
    data_rates_output['version'] = DATA_RATE_OUTPUT_JSON_VER
    data_rates_output['scenario_params'] = sim_case_config['scenario_params']

    viz_output = {}
    viz_output['version'] = VIZ_OUTPUT_JSON_VER
    viz_output['scenario_params'] = sim_case_config['scenario_params']

    # Back to process links copy over
    accesses_data_rates = {}
    viz_data = {}

    include_xlnk_range_in_output = sim_general_config["general_sim_params"]["general_link_params"]['include_xlnk_range_in_output']
    include_dlnk_aer_in_output = sim_general_config["general_sim_params"]["general_link_params"]['include_dlnk_aer_in_output']

    if DATA_RATE_OUTPUT_JSON_VER != "0.3":
        raise NotImplementedError
    if VIZ_OUTPUT_JSON_VER != "0.1":
        raise NotImplementedError

    # Prepare accesses_data_rates
    accesses_data_rates["obs_times"] = rates_output_by_sat_ml["obs"]
    accesses_data_rates["xlnk_times"] = rates_output_by_sat_ml["xlink_update"]
    accesses_data_rates["xlnk_rates"] = rates_output_by_sat_ml["xlink_rates_update"]
    accesses_data_rates["dlnk_times"] = rates_output_by_sat_ml["gslink_update"]
    accesses_data_rates["dlnk_rates"] = rates_output_by_sat_ml["gslink_rates_update"]
    if include_xlnk_range_in_output:
        accesses_data_rates["xlnk_range"] = rates_output_by_sat_ml["xlnkrange_update"]
    if include_dlnk_aer_in_output:
        accesses_data_rates["dlnk_aer"] = rates_output_by_sat_ml["gsaer_update"]
    
    accesses_data_rates["_comment_xlnk_times"]  = "First nesting level is each sat index. Second nesting level is each other sat that can be a cross-link partner. Third nesting level is  each individual cross-link within that sat,xsat pair. For each of the 5-element arrays at that third nesting level,  the first double is the start time of the window (in modified Julian date), the second double is the end time of the window, and the remaining three doubles are intended to be logical values: 1. true if the first satellite is transmitting, 2. true if the second satellite is transmitting, 3. and in the case that they're both transmitting, true if they're transmitting at the same rate over the entire window. Note that this can also be true if only one is transmitting, so it's a necessary but not sufficient condition for symmetry"
    
    accesses_data_rates["_comment_xlnk_rates"]  = "First nesting level is each sat index. Second nesting level is each other sat that can be a cross-link partner (xsat index). Third nesting level is each individual cross-link within that sat,xsat pair. Fourth nesting level is each time point within the cross-link. For each time point there are three fields: 1.  the time (modified Julian date), 2. the transmit data rate in direction sat index -> xsat index, 3. the transmit data rate in direction xsat index -> sat index"

    # Prepare viz_data
    dlnk_partners = viz_output_by_sat_ml["d_part"]
    obs_locations  = viz_output_by_sat_ml["o_locations"]

    viz_data["obs_times_flat"] = viz_output_by_sat_ml["t_o"]
    viz_data["dlnk_times_flat"] = viz_output_by_sat_ml["t_d"]
    viz_data["xlnk_times_flat"] = viz_output_by_sat_ml["t_x"]
    xlnk_partners = viz_output_by_sat_ml["x_part"]
    viz_data["dlnk_rate_history"] = viz_output_by_sat_ml["dlnk_rate_history"]
    viz_data["dlnk_rate_history_epoch"] = \
        sim_case_config['scenario_params']['start_utc'] 
    viz_data["xlnk_rate_history"] = viz_output_by_sat_ml["xlnk_rate_history"]
    viz_data["xlnk_rate_history_epoch"] = \
        sim_case_config['scenario_params']['start_utc'] 
    viz_data["obs_locations"] = obs_locations
    viz_data["dlnk_partners"] = dlnk_partners
    viz_data["xlnk_partners"] = xlnk_partners

    # More from PipelineRunner
    data_rates_output['accesses_data_rates'] = accesses_data_rates
    viz_output['viz_data'] = viz_data

    # We checked in the beginning that the version was 0.3
    data_rates_output['other_data'] = {
        'eclipse_times': orbit_prop_data['accesses_data']['ecl_times']
    }

    return accesses_data_rates, viz_data

def matlab_links_wrapper(sim_case_config, constellation_config,gs_network_config,ops_profile_config,sim_general_config,sat_config,orbit_prop_data,verbose):
    # TODO: There is definitely some redundancy here...should be refactored to match more with Warren's code.
    check_inputs(
        sim_case_config,
        constellation_config,
        gs_network_config,
        ops_profile_config,
        sim_general_config,
        sat_config,
        orbit_prop_data,
        verbose
    )

    accesses_data, params = configure_inputs(
        sim_case_config,
        constellation_config,
        gs_network_config,
        ops_profile_config,
        sim_general_config,
        sat_config,
        orbit_prop_data,
        verbose
    )

    # TODO actually call the link calculator here
    rates_output_by_sat_ml, viz_output_by_sat_ml = links_wrapper.links_wrapper(
        accesses_data, params)

    return configure_outputs(
        rates_output_by_sat_ml,
        viz_output_by_sat_ml,
        sim_case_config,
        constellation_config,
        gs_network_config,
        ops_profile_config,
        sim_general_config,
        sat_config,
        orbit_prop_data,
        verbose
    )
