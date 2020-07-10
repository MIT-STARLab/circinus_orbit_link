#! /usr/bin/env python

##
# Python runner for MATLAB orbit prop pipeline
# @author Kit Kennedy
#

import sys
import time
import os.path
import argparse
import json
from links_wrapper import links_wrapper
import pickle
import numpy as np
import matplotlib.pyplot as plt

#  local repo includes. todo:  make this less hackey
sys.path.append ('..')
try: # First try will work if subrepo circinus_tools is populated, or if prior module imported from elsewhere
    from circinus_tools import io_tools
except ImportError: # Covered by importing orbits above, but don't want to break if that's dropped.
    print("Importing circinus_tools from parent repo...")
    try:
        import sys
        if sys.platform == 'win32':
            sys.path.insert(0, "..\\..\\")
        else:
            sys.path.insert(0, "../../")
        from circinus_tools import io_tools
    except ImportError:
        print("Neither local nor parent-level circinus_tools found.")

REPO_BASE = os.path.abspath(os.pardir)  # os.pardir aka '..'

def assert_lens(arr_1, arr_2, axis):
    assert len(arr_1) == len(arr_2), "axis {} does not match: {} != {}".format(
        axis, len(arr_1), len(arr_2)
    )

# try:
#     import matlab
#     try: # First try will work if subrepo circinus_tools is populated, or if prior module imported from elsewhere
#         from circinus_tools.matlab_if.MatlabIF import MatlabIF
#     except ImportError: # Covered by importing orbits above, but don't want to break if that's dropped.
#         print("Importing circinus_tools from parent repo...")
#         try:
#             import sys
#             sys.path.insert(0, "../../")
#             from circinus_tools.matlab_if.MatlabIF import MatlabIF
#         except ImportError:
#             print("Neither local nor parent-level circinus_tools found.")
#     MATLAB_PIPELINE_ENTRY = os.path.join(REPO_BASE,'matlab_pipeline', 'pipeline_entry')
# except ImportError:
#     print("Matlab module not installed.")



DATA_RATE_OUTPUT_JSON_VER = '0.3' # TODO - vanquish inter-module version checking of this type
VIZ_OUTPUT_JSON_VER = '0.1'

class PipelineRunner:



    def __init__(self,matlab_nestings):
        self.matlabif = None
        self.matlab_nestings =  matlab_nestings

    def process_links(self,sim_case_config,constellation_config,gs_network_config,ops_profile_config,sim_general_config,sat_config,orbit_prop_data,verbose=False):
        """
        :param data: highest-level input json object
        :return: output json with raw orbit prop data
        """

        if verbose:
            print('process_links: convert input to matlab format')

        # matlab-ify the args
        params_ml = {}
        params_ml["scenario_start_utc"] = sim_case_config['scenario_params']['start_utc']
        params_ml["num_sats"] = np.array([constellation_config['constellation_definition']['constellation_params']['num_satellites']])
        params_ml["num_gs"] = matlab.double([gs_network_config['network_definition']['gs_net_params']['num_stations']])
        params_ml["num_targets"] = matlab.double([ops_profile_config['ops_profile_params']['obs_params']['num_targets']])
        sat_id_order= constellation_config['constellation_definition']['constellation_params']['sat_id_order']
        sat_id_prefix= constellation_config['constellation_definition']['constellation_params']['sat_id_prefix']
        gs_id_order= io_tools.make_and_validate_gs_id_order(gs_network_config['network_definition']['gs_net_params'])
        num_sats= constellation_config['constellation_definition']['constellation_params']['num_satellites']

        params_ml["gs_params"] = MatlabIF.deep_convert_python_to_matlab(gs_network_config['network_definition']['gs_net_params'])

        accesses_data_by_sat_ml = {}
        if orbit_prop_data['version'] == "0.3": #only leaving this because it's set in the output of the orbit prop (TODO - drop this versioning too)
            accesses_data_by_sat_ml["obs"]  = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["obs_times"],self.matlab_nestings['obs_times'])
            accesses_data_by_sat_ml["gslink"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["dlnk_times"],self.matlab_nestings['dlnk_times'])
            accesses_data_by_sat_ml["gsaer"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["dlnk_aer"],self.matlab_nestings['dlnk_aer'])
            accesses_data_by_sat_ml["xlink"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["xlnk_times"],self.matlab_nestings['xlnk_times'])
            accesses_data_by_sat_ml["xrange"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["xlnk_range"],self.matlab_nestings['xlnk_range'])
        else:
            raise NotImplementedError

        # if link_inputs['version'] == "0.2":

        params_ml['verbose'] = matlab.logical([sim_general_config["general_sim_params"]["general_link_params"]['matlab_verbose_links']])
        include_xlnk_range_in_output = sim_general_config["general_sim_params"]["general_link_params"]['include_xlnk_range_in_output']
        include_dlnk_aer_in_output = sim_general_config["general_sim_params"]["general_link_params"]['include_dlnk_aer_in_output']

        params_ml["general_link_params"] = MatlabIF.deep_convert_python_to_matlab(sim_general_config["general_sim_params"]["general_link_params"])
        params_ml["lookup_params"] = MatlabIF.deep_convert_python_to_matlab(sim_case_config['scenario_params']["lookup_params"],blind_convert=True)

        # check the first params entry to see if it's using a built in link model

        # TODO - PRI1 - This [0] crap is a hack; get these parameters into the model & move the duplication to the constellation only
        if sat_config['sat_model_definition']['sat_model_params']["sat_link_params"][0]["dlnk_params"]["comm_type"]["built_in"]:
            raise NotImplementedError('Currently the code for built in link models is not included in this repository')
        if sat_config['sat_model_definition']['sat_model_params']["sat_link_params"][0]["xlnk_params"]["comm_type"]["built_in"]:
            raise NotImplementedError('Currently the code for built in link models is not included in this repository')

        slp_holder = {
            "sat_id_prefix": constellation_config['constellation_definition']['constellation_params']['sat_id_prefix'],
            "sat_ids": constellation_config['constellation_definition']['constellation_params']['sat_ids'],
            **sat_config['sat_model_definition']['sat_model_params']['sat_link_params'][0]
        }

        sat_config['sat_model_definition']['sat_model_params']['sat_link_params'] = [slp_holder]

        # in the case that this is default, then we need to grab a list of all the satellite IDs. We'll take this from all of the satellite IDs found in the orbit parameters
        if sat_id_order == 'default': # TODO - PRI1 - This matlab breakdown is dependent on looking up by the zhou name in sat_link_params>dlink_params>comm_type>name. Make it not so.
            sat_link_params, all_sat_ids = io_tools.unpack_sat_entry_list(  sat_config['sat_model_definition']['sat_model_params']["sat_link_params"],force_duplicate =  True)

        #  make the satellite ID order. if the input ID order is default, then will assume that the order is the same as all of the IDs passed as argument
        sat_id_order = io_tools.make_and_validate_sat_id_order(sat_id_order,sat_id_prefix,num_sats,all_sat_ids)
        io_tools.validate_ids(validator=sat_id_order,validatee=all_sat_ids)

        sat_link_params_sorted = io_tools.sort_input_params_by_sat_IDs(sat_link_params,sat_id_order)

        params_ml["sat_link_params"] = MatlabIF.deep_convert_python_to_matlab(sat_link_params_sorted)

        def get_sat_indx_str(sat_ID):
            return str(sat_id_order.index(str(sat_ID)))
        def get_gs_indx_str(gs_ID):
            return str(gs_id_order.index(str(gs_ID)))
        def convert_list_elems(lst,f_conv):
            return [f_conv(elem) for elem in lst]

        dlnk_direc_disabled_gs_ID_by_sat_IDstr = ops_profile_config['ops_profile_params']["link_disables"]['dlnk_direc_disabled_gs_ID_by_sat_IDstr']
        xlnk_direc_disabled_xsat_ID_by_sat_IDstr = ops_profile_config['ops_profile_params']["link_disables"]['xlnk_direc_disabled_xsat_ID_by_sat_IDstr']
        #  convert all elements over to strings so they can be passed to Matlab.  need to prepend a string to the satellite IDs so that they can be valid fields in a Matlab struct
        dlnk_direc_disabled_gs_indx_by_sat_indx = {'indx'+get_sat_indx_str(sat_ID): convert_list_elems(disable_IDs,get_gs_indx_str) for sat_ID,disable_IDs in dlnk_direc_disabled_gs_ID_by_sat_IDstr.items()}
        xlnk_direc_disabled_xsat_indx_by_sat_indx = {'indx'+get_sat_indx_str(sat_ID): convert_list_elems(disable_IDs,get_sat_indx_str) for sat_ID,disable_IDs in xlnk_direc_disabled_xsat_ID_by_sat_IDstr.items()}

        params_ml['link_disables'] = {
            "dlnk_direc_disabled_gs_indx_by_sat_indx": dlnk_direc_disabled_gs_indx_by_sat_indx,
            "xlnk_direc_disabled_xsat_indx_by_sat_indx": xlnk_direc_disabled_xsat_indx_by_sat_indx,
        }


        if not self.matlabif:
            self.matlabif = MatlabIF(matlab_ver=sim_general_config["general_sim_params"]["matlab_version"],paths=[MATLAB_PIPELINE_ENTRY])

        if verbose:
            print('process_links: call matlab')

        (rates_output_by_sat_ml,viz_output_by_sat_ml) = links_wrapper.links_wrapper(accesses_data_by_sat_ml, params_ml)

        accesses_data_rates = {}
        viz_data = {}

        if verbose:
            print('process_links: convert outputs')

        if DATA_RATE_OUTPUT_JSON_VER == "0.3":
            accesses_data_rates["obs_times"]  = MatlabIF.deep_convert_matlab_to_python(rates_output_by_sat_ml["obs"])
            accesses_data_rates["xlnk_times"]  = MatlabIF.deep_convert_matlab_to_python(rates_output_by_sat_ml["xlink_update"])
            accesses_data_rates["xlnk_rates"] = MatlabIF.deep_convert_matlab_to_python(rates_output_by_sat_ml["xlink_rates_update"])
            accesses_data_rates["dlnk_times"] = MatlabIF.deep_convert_matlab_to_python(rates_output_by_sat_ml["gslink_update"])
            accesses_data_rates["dlnk_rates"] = MatlabIF.deep_convert_matlab_to_python(rates_output_by_sat_ml["gslink_rates_update"])
            if include_xlnk_range_in_output:
                accesses_data_rates["xlnk_range"] = MatlabIF.deep_convert_matlab_to_python(rates_output_by_sat_ml["xlnkrange_update"])
            if include_dlnk_aer_in_output:
                accesses_data_rates["dlnk_aer"] = MatlabIF.deep_convert_matlab_to_python(rates_output_by_sat_ml["gsaer_update"])

            accesses_data_rates["_comment_xlnk_times"]  = "First nesting level is each sat index. Second nesting level is each other sat that can be a cross-link partner. Third nesting level is  each individual cross-link within that sat,xsat pair. For each of the 5-element arrays at that third nesting level,  the first double is the start time of the window (in modified Julian date), the second double is the end time of the window, and the remaining three doubles are intended to be logical values: 1. true if the first satellite is transmitting, 2. true if the second satellite is transmitting, 3. and in the case that they're both transmitting, true if they're transmitting at the same rate over the entire window. Note that this can also be true if only one is transmitting, so it's a necessary but not sufficient condition for symmetry"

            accesses_data_rates["_comment_xlnk_rates"]  = "First nesting level is each sat index. Second nesting level is each other sat that can be a cross-link partner (xsat index). Third nesting level is each individual cross-link within that sat,xsat pair. Fourth nesting level is each time point within the cross-link. For each time point there are three fields: 1.  the time (modified Julian date), 2. the transmit data rate in direction sat index -> xsat index, 3. the transmit data rate in direction xsat index -> sat index"
        else:
            raise NotImplementedError

        if VIZ_OUTPUT_JSON_VER == "0.1": # TODO - this goes straight to another module, so excusing the version numbering but also needs to go
            viz_data["obs_times_flat"]  = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["t_o"],self.matlab_nestings['obs_times_flat'])
            obs_locations  = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["o_locations"],self.matlab_nestings['obs_locations'])
            viz_data["dlnk_times_flat"] = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["t_d"],self.matlab_nestings['dlnk_times_flat'])
            dlnk_partners = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["d_part"],self.matlab_nestings['dlnk_partners'])
            viz_data["xlnk_times_flat"] = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["t_x"],self.matlab_nestings['xlnk_times_flat'])
            xlnk_partners = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["x_part"],self.matlab_nestings['xlnk_partners'])
            viz_data["dlnk_rate_history"] = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["dlnk_rate_history"])
            viz_data["dlnk_rate_history_epoch"] = sim_case_config['scenario_params']['start_utc']
            viz_data["xlnk_rate_history"] = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["xlnk_rate_history"])
            viz_data["xlnk_rate_history_epoch"] = sim_case_config['scenario_params']['start_utc']

            viz_data["obs_locations"] = MatlabIF.convert_matlab_indexing_to_python (obs_locations)
            viz_data["dlnk_partners"] = MatlabIF.convert_matlab_indexing_to_python (dlnk_partners)
            viz_data["xlnk_partners"] = MatlabIF.convert_matlab_indexing_to_python (xlnk_partners)

        else:
            raise NotImplementedError



        return accesses_data_rates, viz_data


    def get_other_data_rates_output( self,orbit_prop_data):

        other_data = {}

        if orbit_prop_data['version'] == "0.3":
            other_data['eclipse_times'] = orbit_prop_data['accesses_data']['ecl_times']
        else:
            raise NotImplementedError

        return other_data

    def run(self, data):
        """
        Run orbit propagation pipeline element using the inputs supplied per input.json schema. Formats the high level output json and calls various subcomponents for processing

        :param data: input json per input.json schema
        :return: output json per output.json schema
        """

        orbit_prop_data = data['orbit_prop_data']

        sim_case_config = data['sim_case_config'] # TODO; why pack it into data and then re-unpack it once in run?
        constellation_config = data['constellation_config']
        gs_network_config = data['gs_network_config']
        ops_profile_config = data['ops_profile_config']
        sim_general_config = data['sim_general_config']
        sat_config = data['sat_config']


        # define orbit prop outputs json
        data_rates_output = {}
        data_rates_output['version'] = DATA_RATE_OUTPUT_JSON_VER
        data_rates_output['scenario_params'] = sim_case_config['scenario_params']

        viz_output = {}
        viz_output['version'] = VIZ_OUTPUT_JSON_VER
        viz_output['scenario_params'] = sim_case_config['scenario_params']




        out_stuff = self.process_links(
            sim_case_config,
            constellation_config,
            gs_network_config,
            ops_profile_config,
            sim_general_config,
            sat_config,
            orbit_prop_data,
            verbose=True)

        data_rates_output['accesses_data_rates'] = out_stuff[0]
        data_rates_output['other_data'] = self.get_other_data_rates_output (orbit_prop_data)
        viz_output['viz_data'] = out_stuff[1]

        return data_rates_output, viz_output



def main():
    ap = argparse.ArgumentParser(description='orbit communication link calculations')

    # Using default arguments for providing the default location that would be used if not w/in larger CIRCINUS repo.
    # TODO - When we settle on our nominal models, will make this example directory.
    ap.add_argument('--inputs_location',
                    type=str,
                    default=os.path.join(REPO_BASE,'example_input'),
                    help='specify directory in which to find input and config files')

    ap.add_argument('--case_name',
                    type=str,
                    default=os.path.join(REPO_BASE,'nom_case'),
                    help='specify name of case to be used for calculations')

    # TODO: Add default example input w/in this repo

    ap.add_argument('--prop_data_inputs',
                    type=str,
                    default='orbit_prop_data.json',
                    help='specify orbit link inputs file')

    ap.add_argument('--link_calculator_language',
                    type=str,
                    default='python',
                    help='specify language to use, either: matlab or python available')

    args = ap.parse_args()

    # ------- Filenames ------ #
    sim_case_config_FILENAME = args.inputs_location+'/cases/'+args.case_name+'/sim_case_config.json'
    constellation_config_FILENAME = args.inputs_location+'/cases/'+args.case_name+'/constellation_config.json'
    gs_network_config_FILENAME = args.inputs_location+'/cases/'+args.case_name+'/ground_station_network_config.json'
    ops_profile_config_FILENAME = args.inputs_location+'/cases/'+args.case_name+'/operational_profile_config.json'
    sim_general_config_FILENAME = args.inputs_location+'/general_config/sim_general_config.json'
    orbit_prop_data_FILENAME = args.inputs_location+'/cases/'+args.case_name+'/autogen_files/orbit_prop_data.json'

    # Winsows users check
    import sys
    if sys.platform == 'win32':
        sim_case_config_FILENAME.replace('/','\\')
        sim_general_config_FILENAME.replace('/','\\')
        constellation_config_FILENAME.replace('/','\\')
        gs_network_config_FILENAME.replace('/','\\')
        ops_profile_config_FILENAME.replace('/','\\')


    # ------- Check for file validity ------ #
    sim_case_config_EXISTS = os.path.isfile(sim_case_config_FILENAME)
    constellation_config_EXISTS = os.path.isfile(constellation_config_FILENAME)
    gs_network_config_EXISTS = os.path.isfile(gs_network_config_FILENAME)
    ops_profile_config_EXISTS = os.path.isfile(ops_profile_config_FILENAME)
    sim_general_config_EXISTS = os.path.isfile(sim_general_config_FILENAME)
    orbit_prop_data_EXISTS = os.path.isfile(orbit_prop_data_FILENAME)
    sat_config_EXISTS = False # Covered below

    # ------- Extract another filename to check from one of the files above ------- #
    # Constellation Config points us to the satellite models being used - # Todo: handle multiple sat models
    if(constellation_config_EXISTS):
        with open(constellation_config_FILENAME,'r') as f:
            constellation_config = json.load(f)
        sat_ref_model_name = constellation_config["constellation_definition"]["default_sat_ref_model_name"] # Only handling default at the moment
        sat_config_FILENAME = args.inputs_location+'/reference_model_definitions/sat_refs/'+sat_ref_model_name+'.json'
        sat_config_EXISTS = os.path.isfile(sat_config_FILENAME)

    if(not sim_case_config_EXISTS
        or not constellation_config_EXISTS
        or not gs_network_config_EXISTS
        or not ops_profile_config_EXISTS
        or not sim_general_config_EXISTS
        or not sat_config_EXISTS
        or not orbit_prop_data_EXISTS):

        print("Some required configuration files may not provided:")
        print("  in inputs/cases/{}:".format(args.case_name))
        print("\tconstellation_config.json              {}".format(("[missing]" if not constellation_config_EXISTS else "")))
        print("\tground_station_network_config.json     {}".format(("[missing]" if not gs_network_config_EXISTS else "")))
        print("\toperational_profile_config.json        {}".format(("[missing]" if not ops_profile_config_EXISTS else "")))
        print("\tsim_case_config.json                   {}".format(("[missing]" if not sim_case_config_EXISTS else "")))

        print("  in inputs/general_config:")
        print("\tsim_general_config.json                {}".format(("[missing]" if not sim_general_config_EXISTS else "")))

        # TODO - this should be abstracted out via the constellation file
        print("  in inputs/reference_model_definitions/sat_refs:")
        print("\tnominal_original_sat.json              {}".format(("[missing]" if not sat_config_EXISTS else "")))

        print("  in inputs/cases/{}/autogen_files:".format(args.case_name))
        print("\torbit_prop_data.json                   {}".format(("[missing] (created by running the orbit_propagation module prior to this.)" if not orbit_prop_data_EXISTS else "")))

        print("Provide any missing input files. Terminating.\n")
        exit(1)



    # -------- CASE SPECIFIC CONFIGURATION INPUTS -------- #
    with open(sim_case_config_FILENAME,'r') as f:
        sim_case_config = json.load(f)

    # with open(constellation_config_FILENAME,'r') as f:    # Covered by check above where we need to extract the sat model name

    with open(gs_network_config_FILENAME,'r') as f:
        gs_network_config = json.load(f)

    with open(ops_profile_config_FILENAME,'r') as f:
        ops_profile_config = json.load(f)


    # -------- GENERAL CONFIGURATION INPUTS -------- #
    with open(sim_general_config_FILENAME,'r') as f:
        sim_general_config = json.load(f)

    with open(sat_config_FILENAME,'r') as f:
        sat_config = json.load(f)



    # -------- AUTO-GENERATED INPUTS (from prior modules) -------- #
    with open(orbit_prop_data_FILENAME,'r') as f:
        orbit_prop_data = json.load(f)


    data = {
        "sim_case_config":sim_case_config,
        "constellation_config":constellation_config,
        "gs_network_config":gs_network_config,
        "ops_profile_config":ops_profile_config,

        "sim_general_config":sim_general_config,
        "sat_config":sat_config,

        # Auto-gen:
        "orbit_prop_data": orbit_prop_data,
    }

    a = time.time()
    if args.link_calculator_language == 'matlab':
        print("Language is matlab!")
        MATLAB_SINGLE_NESTED_MAT = [matlab.double]
        MATLAB_DOUBLE_NESTED_MAT = [[matlab.double]]
        MATLAB_TRIPLE_NESTED_MAT = [[[matlab.double]]]
        MATLAB_QUADRUPLE_NESTED_MAT = [[[[matlab.double]]]]
        MATLAB_DOUBLE_NESTED_FLAT = [['flat_list']]
        MATLAB_DOUBLE_NESTED_VAL = [['value']]

        matlab_nesting = {
                'obs_times': MATLAB_DOUBLE_NESTED_MAT,
                'dlnk_times': MATLAB_DOUBLE_NESTED_MAT,
                'xlnk_times': MATLAB_DOUBLE_NESTED_MAT,
                'obs_aer': MATLAB_TRIPLE_NESTED_MAT,
                'xlnk_range': MATLAB_TRIPLE_NESTED_MAT,
                'dlnk_aer': MATLAB_TRIPLE_NESTED_MAT,
                'ecl_times': MATLAB_SINGLE_NESTED_MAT,
                'obs_times_flat': MATLAB_DOUBLE_NESTED_FLAT,
                'obs_locations': MATLAB_DOUBLE_NESTED_VAL,
                'dlnk_times_flat': MATLAB_DOUBLE_NESTED_FLAT,
                'dlnk_partners': MATLAB_DOUBLE_NESTED_VAL,
                'xlnk_times_flat': MATLAB_DOUBLE_NESTED_FLAT,
                'xlnk_partners': MATLAB_DOUBLE_NESTED_VAL,
            }
        pr = PipelineRunner(matlab_nesting)
        output = pr.run(data)
    elif args.link_calculator_language == 'python':
        from simple_link_calc import py_links_wrapper
        output = py_links_wrapper(data)
        # # accesses_data_rates, viz_data = py_links_wrapper(data)

        # out = []
        # with open("/home/dani/Documents/star_lab/orig_io/runner_orbitlink_saveoff/"
        #         "process_links_output.pkl", "rb") as f:
        #     out = pickle.load(f)

        # real_accesses_data_rates = out[0]
        # real_viz_data = out[1]
        # assert real_accesses_data_rates['_comment_xlnk_rates'] == accesses_data_rates['_comment_xlnk_rates']


        # assert real_accesses_data_rates['_comment_xlnk_times'] == accesses_data_rates['_comment_xlnk_times']

        # def dump_dlnk_rates(lnk_rates, name):
        #     print("------------ " + name.upper() + "RATES ------------")
        #     print(
        #         "len axis 0 (num sats)  : {}\n"
        #         "len axis 1 (num gs)    : {}\n"
        #         "len axis 2 (num lnks) : {}\n"
        #         "len axis 3: {}\n"
        #         "len axis 4: {}".format(
        #         len(lnk_rates), len(lnk_rates[-1]), len(lnk_rates[-1][-1]), len(lnk_rates[-1][-1][-1]),
        #         len(lnk_rates[-1][-1][-1][-1]),
        #     ))
        #     print(name + "[0][0][0][0]:", lnk_rates[0][0][0][0])
        #     print(name + "[-2][-2][-2][-2]:", lnk_rates[-2][-2][-2][-2])
        #     print(name + "[-1][-1][-1][-1]:", lnk_rates[-1][-1][-1][-1])

        # try:
        #     assert real_accesses_data_rates['dlnk_rates'] == accesses_data_rates['dlnk_rates']
        # except AssertionError:
        #     print("=========== DOWNLINK RATES ============")
        #     dump_dlnk_rates(real_accesses_data_rates['dlnk_rates'], "real dlnk")
        #     dump_dlnk_rates(accesses_data_rates['dlnk_rates'], "mine dlnk")
        #     raise

        # print("Dlnk rates match!")

        # def dump_times(lnk_times, name):
        #     print("------------ " + name.upper() + " TIMES ------------")
        #     print(
        #         "len axis 0 (num sats)  : {}\n"
        #         "len axis 1 (num gs)    : {}\n"
        #         "len axis 2 (num lnks) : {}\n"
        #         "len axis 3: {}".format(
        #         len(lnk_times), len(lnk_times[-1]), len(lnk_times[-1][-1]), len(lnk_times[-1][-1][-1]),
        #     ))
        #     print(name + "[0][0][0][0]:", lnk_times[0][0][0][0])
        #     print(name + "[-2][-2][-2][-2]:", lnk_times[-2][-2][-2][-2])
        #     print(name + "[-1][-1][-1][-1]:", lnk_times[-1][-1][-1][-1])

        # def compare_dlnk_times(real, mine):
        #     assert len(real) == len(mine), \
        #         "len(real) = {} != len(mine) = {}".format(len(real), len(mine))
        #     for sat_num in range(len(real)):
        #         assert len(real[sat_num]) == len(mine[sat_num]), \
        #             "len(real[{}]) = {} != len(mine[{}]) = {}".format(sat_num,
        #             len(real), sat_num, len(mine))
        #         for gs_num in range(len(real[sat_num])):
        #             assert len(real[sat_num]) == len(mine[sat_num]), \
        #                 "len(real[{}][{}]) = {} != len(mine[{}][{}]) = {}".format(
        #                     sat_num, gs_num, len(real), sat_num, gs_num, len(mine))
        #             for dlnk_num in range(len(real[sat_num][gs_num])):
        #                 if real[sat_num][gs_num][dlnk_num] != mine[sat_num][gs_num][dlnk_num]:
        #                     print("real[{}][{}][{}] = {} != mine[{}][{}][{}] = {}".format(
        #                         sat_num, gs_num, dlnk_num, real[sat_num][gs_num][dlnk_num], sat_num, gs_num, dlnk_num, mine[sat_num][gs_num][dlnk_num]
        #                     ))
        #     print("All lengths match!")

        # try:
        #     assert real_accesses_data_rates['dlnk_times'] == accesses_data_rates['dlnk_times']
        # except AssertionError:
        #     print("=========== DOWNLINK TIMES ============")
        #     dump_times(real_accesses_data_rates['dlnk_times'], "real")
        #     dump_times(accesses_data_rates['dlnk_times'], "mine")
        #     compare_dlnk_times(real_accesses_data_rates['dlnk_times'], accesses_data_rates['dlnk_times'])
        #     raise

        # print("Dlnk times match!")

        # assert real_accesses_data_rates['obs_times'] == accesses_data_rates['obs_times']

        # print("Obs times match!")

        # def dump_xlnk_rates(lnk_rates, name):
        #     non_zero_indxs = []
        #     print("------------ " + name.upper() + " RATES ------------")
        #     print(
        #         "len axis 0: {}\n"
        #         "len axis 1: {}".format(
        #         len(lnk_rates), len(lnk_rates[-1])
        #     ))
        #     print(name + "[0][0]:", lnk_rates[0][0])
        #     for i, sup_list in enumerate(lnk_rates):
        #         for j, l in enumerate(sup_list):
        #             if len(l) == 0:
        #                 continue
        #             non_zero_indxs.append((i, j))
        #             print("len(lnk_rates[{}][{}]: {}".format(i, j, len(l)))
        #     return non_zero_indxs

        # def compare_xlnk_rates(real, mine, non_zero_indxs):
        #     for ind_0, ind_1 in non_zero_indxs:
        #         for i, (r_i_arr, m_i_arr) in enumerate(zip(real[ind_0][ind_1],
        #                 mine[ind_0][ind_1])):
        #             assert_lens(r_i_arr, m_i_arr, i)
        #             for j, (r_j_arr, m_j_arr) in enumerate(zip(r_i_arr, m_i_arr)):
        #                 assert_lens(r_j_arr, m_j_arr, j)
        #                 for k, (r_k, m_k) in enumerate(zip(r_j_arr, m_j_arr)):
        #                     if r_k != m_k:
        #                         print(
        #                             "real[{}][{}][{}][{}][{}] = {} != "
        #                             "mine[{}][{}][{}][{}][{}] = {}".format(
        #                                 ind_0, ind_1, i, j, k, r_k,
        #                                 ind_0, ind_1, i, j, k, m_k
        #                             )
        #                         )
        #                     else:
        #                         print("Equal at [{}][{}][{}]".format(i, j, k))
        #     print("All dims match!")

        # try:
        #     assert real_accesses_data_rates['xlnk_rates'] == accesses_data_rates['xlnk_rates']
        # except AssertionError:
        #     print("============ CROSSLINK RATES ============")
        #     rnz = dump_xlnk_rates(real_accesses_data_rates['xlnk_rates'], "real xlnk")
        #     mnz = dump_xlnk_rates(accesses_data_rates['xlnk_rates'], "mine xlnk")
        #     assert rnz == mnz
        #     compare_xlnk_rates(real_accesses_data_rates['xlnk_rates'],
        #         accesses_data_rates['xlnk_rates'], rnz)
        #     raise

        # print("Xlnk rates match!")

        # assert real_accesses_data_rates['xlnk_times'] == accesses_data_rates['xlnk_times']

        # print("Xlnk times match!")

        # assert real_accesses_data_rates == accesses_data_rates
    
        # print("*** All access data match! ***")

        # for key in sorted(real_viz_data.keys()):
        #     if key == 'dlnk_partners' or key == 'dlnk_times_flat':
        #         assert len(real_viz_data[key]) == \
        #             len(viz_data[key]), \
        #                 "real len ({}) != my len ({})".format(
        #                     len(real_viz_data[key]),
        #                     len(viz_data[key])
        #                 )
        #         for i, (r_val, m_val) in enumerate(zip(real_viz_data[key],
        #                 viz_data[key])):
        #             assert len(r_val) == len(m_val), '|r[{}]| = {} != |m[{}]| = {}'.format(
        #                 i, len(r_val), i, len(m_val)
        #             )
        #             for j, (rr_val, mm_val) in enumerate(zip(r_val, m_val)):
        #                 if rr_val != mm_val:
        #                     print("real[{}][{}] = {} != mine[{}][{}] = {}".format(
        #                         i, j, rr_val, i, j, mm_val
        #                     ))
        #     assert real_viz_data[key] == viz_data[key],  \
        #         key + " in viz data do not match!"
        #     print(key, "in viz data match!")

        # print("Calc done, outputs equal!")
    else:
        raise NotImplementedError('only python and matlab are available language options')

    b = time.time()


    with open(args.inputs_location+'/cases/'+args.case_name+'/autogen_files/data_rates_output.json','w') as f:
        json.dump(output[0], f, indent=4, separators=(',', ': '))
    with open(args.inputs_location+'/cases/'+args.case_name+'/autogen_files/sat_link_history.json','w') as f:
        json.dump(output[1], f, indent=4, separators=(',', ': '))

    print('run time: %f'%(b-a))

if __name__ == "__main__":
    main()
