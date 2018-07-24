#! /usr/bin/env python

##
# Python runner for MATLAB orbit prop pipeline
# @author Kit Kennedy
# 

import sys
import time
import os.path
import matlab
import argparse

#  local repo includes. todo:  make this less hackey
sys.path.append ('..')
from circinus_tools import io_tools
from circinus_tools.matlab_if import MatlabIF

REPO_BASE = os.path.abspath(os.pardir)  # os.pardir aka '..'
MATLAB_PIPELINE_ENTRY = os.path.join(REPO_BASE,'matlab_pipeline', 'pipeline_entry')

DATA_RATE_OUTPUT_JSON_VER = '0.3'
VIZ_OUTPUT_JSON_VER = '0.1'

MATLAB_SINGLE_NESTED_MAT = [matlab.double]
MATLAB_DOUBLE_NESTED_MAT = [[matlab.double]]
MATLAB_TRIPLE_NESTED_MAT = [[[matlab.double]]]
MATLAB_QUADRUPLE_NESTED_MAT = [[[[matlab.double]]]]
MATLAB_DOUBLE_NESTED_FLAT = [['flat_list']]
MATLAB_DOUBLE_NESTED_VAL = [['value']]

class PipelineRunner:

    matlab_nestings =  {
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

    def __init__(self):
        self.matlabif = None

    def process_links(self,orbit_prop_inputs,orbit_prop_data,link_inputs,verbose=False):
        """

        :param data: highest-level input json object
        :return: output json with raw orbit prop data
        """

        if verbose:
            print('process_links: convert input to matlab format')

        # matlab-ify the args
        params_ml = {}
        if orbit_prop_inputs['version'] == "0.8":
            params_ml["scenario_start_utc"] = orbit_prop_inputs['scenario_params']['start_utc']
            params_ml["num_sats"] = matlab.double([orbit_prop_inputs['sat_params']['num_satellites']])
            params_ml["num_gs"] = matlab.double([orbit_prop_inputs['gs_params']['num_stations']])
            params_ml["num_targets"] = matlab.double([orbit_prop_inputs['obs_params']['num_targets']])
            sat_id_order= orbit_prop_inputs['sat_params']['sat_id_order']
            sat_id_prefix= orbit_prop_inputs['sat_params']['sat_id_prefix']
            gs_id_order= io_tools.make_and_validate_gs_id_order(orbit_prop_inputs['gs_params'])
            num_sats= orbit_prop_inputs['sat_params']['num_satellites']

            params_ml["gs_params"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_inputs['gs_params'])
        else:
            raise NotImplementedError

        accesses_data_by_sat_ml = {}
        if orbit_prop_data['version'] == "0.3":
            accesses_data_by_sat_ml["obs"]  = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["obs_times"],self.matlab_nestings['obs_times'])
            accesses_data_by_sat_ml["gslink"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["dlnk_times"],self.matlab_nestings['dlnk_times'])
            accesses_data_by_sat_ml["gsaer"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["dlnk_aer"],self.matlab_nestings['dlnk_aer'])
            accesses_data_by_sat_ml["xlink"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["xlnk_times"],self.matlab_nestings['xlnk_times'])
            accesses_data_by_sat_ml["xrange"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["xlnk_range"],self.matlab_nestings['xlnk_range'])
            # accesses_data_by_sat_ml["obs_aer"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["obs_aer"])
            # accesses_data_by_sat_ml["ecl_times"] = MatlabIF.deep_convert_python_to_matlab(orbit_prop_data["accesses_data"]["ecl_times"])
        else:
            raise NotImplementedError

        if link_inputs['version'] == "0.2":
            
            params_ml['verbose'] = matlab.logical([link_inputs["general_link_params"]['matlab_verbose_links']])
            include_xlnk_range_in_output = link_inputs["general_link_params"]['include_xlnk_range_in_output']
            include_dlnk_aer_in_output = link_inputs["general_link_params"]['include_dlnk_aer_in_output']

            params_ml["general_link_params"] = MatlabIF.deep_convert_python_to_matlab(link_inputs["general_link_params"])
            params_ml["lookup_params"] = MatlabIF.deep_convert_python_to_matlab(link_inputs["lookup_params"],blind_convert=True)

            # check the first params entry to see if it's using a built in link model
            if link_inputs["sat_link_params"][0]["dlnk_params"]["comm_type"]["built_in"]:
                raise NotImplementedError('Currently the code for built in link models is not included in this repository')
            if link_inputs["sat_link_params"][0]["xlnk_params"]["comm_type"]["built_in"]:
                raise NotImplementedError('Currently the code for built in link models is not included in this repository')

            # in the case that this is default, then we need to grab a list of all the satellite IDs. We'll take this from all of the satellite IDs found in the orbit parameters
            if sat_id_order == 'default':
                sat_link_params, all_sat_ids = io_tools.unpack_sat_entry_list(  link_inputs["sat_link_params"],force_duplicate =  True)
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

            dlnk_direc_disabled_gs_ID_by_sat_IDstr = link_inputs["link_disables"]['dlnk_direc_disabled_gs_ID_by_sat_IDstr']
            xlnk_direc_disabled_xsat_ID_by_sat_IDstr = link_inputs["link_disables"]['xlnk_direc_disabled_xsat_ID_by_sat_IDstr']
            #  convert all elements over to strings so they can be passed to Matlab.  need to prepend a string to the satellite IDs so that they can be valid fields in a Matlab struct
            dlnk_direc_disabled_gs_indx_by_sat_indx = {'indx'+get_sat_indx_str(sat_ID): convert_list_elems(disable_IDs,get_gs_indx_str) for sat_ID,disable_IDs in dlnk_direc_disabled_gs_ID_by_sat_IDstr.items()}
            xlnk_direc_disabled_xsat_indx_by_sat_indx = {'indx'+get_sat_indx_str(sat_ID): convert_list_elems(disable_IDs,get_sat_indx_str) for sat_ID,disable_IDs in xlnk_direc_disabled_xsat_ID_by_sat_IDstr.items()}

            params_ml['link_disables'] = {
                "dlnk_direc_disabled_gs_indx_by_sat_indx": dlnk_direc_disabled_gs_indx_by_sat_indx,
                "xlnk_direc_disabled_xsat_indx_by_sat_indx": xlnk_direc_disabled_xsat_indx_by_sat_indx,
            }

            # from circinus_tools import debug_tools
            # debug_tools.debug_breakpt()

        else:
            raise NotImplementedError



        if not self.matlabif:
            self.matlabif = MatlabIF(paths=[MATLAB_PIPELINE_ENTRY])
            
        if verbose:
            print('process_links: call matlab')

        (rates_output_by_sat_ml,viz_output_by_sat_ml) = self.matlabif.call_mfunc(
                'links_wrapper',
                accesses_data_by_sat_ml,
                params_ml,
                nargout=2)

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

        if VIZ_OUTPUT_JSON_VER == "0.1":
            viz_data["obs_times_flat"]  = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["t_o"],self.matlab_nestings['obs_times_flat'])
            obs_locations  = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["o_locations"],self.matlab_nestings['obs_locations'])
            viz_data["dlnk_times_flat"] = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["t_d"],self.matlab_nestings['dlnk_times_flat'])
            dlnk_partners = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["d_part"],self.matlab_nestings['dlnk_partners'])
            viz_data["xlnk_times_flat"] = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["t_x"],self.matlab_nestings['xlnk_times_flat'])
            xlnk_partners = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["x_part"],self.matlab_nestings['xlnk_partners'])
            viz_data["dlnk_rate_history"] = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["dlnk_rate_history"])
            viz_data["dlnk_rate_history_epoch"] = orbit_prop_inputs['scenario_params']['start_utc']
            viz_data["xlnk_rate_history"] = MatlabIF.deep_convert_matlab_to_python(viz_output_by_sat_ml["xlnk_rate_history"])
            viz_data["xlnk_rate_history_epoch"] = orbit_prop_inputs['scenario_params']['start_utc']

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

    def run(self,data):
        """
        Run orbit propagation pipeline element using the inputs supplied per input.json schema. Formats the high level output json and calls various subcomponents for processing

        :param data: input json per input.json schema
        :return: output json per output.json schema
        """

        orbit_prop_data = data['orbit_prop_data']
        orbit_prop_inputs = data['orbit_prop_inputs']

        # define orbit prop outputs json
        data_rates_output = {}
        data_rates_output['version'] = DATA_RATE_OUTPUT_JSON_VER
        data_rates_output['scenario_params'] = orbit_prop_inputs['scenario_params']

        viz_output = {}
        viz_output['version'] = VIZ_OUTPUT_JSON_VER
        viz_output['scenario_params'] = orbit_prop_inputs['scenario_params']

        link_inputs=data['link_inputs']
        out_stuff =self.process_links(orbit_prop_inputs,orbit_prop_data,link_inputs,verbose=True)

        data_rates_output['accesses_data_rates'] = out_stuff[0]
        data_rates_output['other_data'] = self.get_other_data_rates_output (orbit_prop_data)
        viz_output['viz_data'] = out_stuff[1]

        return data_rates_output, viz_output



if __name__ == "__main__":

    ap = argparse.ArgumentParser(description='orbit communication link calculations')
    ap.add_argument('--prop_inputs_file',
                    type=str,
                    default='orbit_prop_inputs.json',
                    help='specify orbit propagation inputs file')

    ap.add_argument('--link_inputs_file',
                    type=str,
                    default='orbit_link_inputs.json',
                    help='specify orbit link inputs file')

    args = ap.parse_args()

    pr = PipelineRunner()

    import json

    with open(os.path.join(REPO_BASE,args.prop_inputs_file),'r') as f:
        orbit_prop_inputs = json.load(f)

    with open(os.path.join(REPO_BASE,'crux/config/examples/orbit_prop_data.json'),'r') as f:
    # with open(os.path.join(REPO_BASE,'crux/config/examples/orbit_prop_data_ex.json'),'r') as f:
        orbit_prop_data = json.load(f)

    with open(os.path.join(REPO_BASE,args.link_inputs_file),'r') as f:
        link_inputs = json.load(f)

    data = {
        "orbit_prop_data": orbit_prop_data,
        "orbit_prop_inputs": orbit_prop_inputs,
        "link_inputs":  link_inputs
    }

    a = time.time()
    output = pr.run(data)
    b = time.time()

    with open('data_rates_output.json','w') as f:
        json.dump(output[0],f)
    with open('sat_link_history.json','w') as f:
        json.dump(output[1],f)

    print('run time: %f'%(b-a))
