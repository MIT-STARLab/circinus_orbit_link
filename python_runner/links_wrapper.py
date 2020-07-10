import links_calculator
import os
import numpy as np


def reshape_inputs(in_by_sat, dim_2_size):
    """
    Copied from circinus_orbit_link_public/matlab_pipeline/matlab_tools/reshape_inputs.m
    """
    return [i[:dim_2_size] for i in in_by_sat]


def links_wrapper(accesses_data_by_sat, params):
    if params["verbose"]:  # TODO: Fix this
        print("running links wrapper")

    # Find bas directory for matlab pipeline code
    # file_dir = os.path.dirname(os.path.realpath(__file__))  # Current dir
    # base_directory = os.path.join(file_dir, "..")
    # base_directory = os.path.join(base_directory, "matlab_tools")
    # TODO: Make sure this is the right directory

    accesses_data = {}
    accesses_data["obs"] = reshape_inputs(accesses_data_by_sat["obs"],
                                          params["num_targets"])
    accesses_data["gslink"] = reshape_inputs(accesses_data_by_sat["gslink"],
                                             params["num_gs"])
    accesses_data["gsaer"] = reshape_inputs(accesses_data_by_sat["gsaer"],
                                            params["num_gs"])
    accesses_data["xlink"] = reshape_inputs(accesses_data_by_sat["xlink"],
                                            params["num_sats"])
    accesses_data["xrange"] = reshape_inputs(accesses_data_by_sat["xrange"],
                                             params["num_sats"])

    print("Links wrapper called!")

    ret = links_calculator.links_calculator(
        accesses_data, params, params['verbose'])

    print("Links calculator didn't crash!")

    rates_output, viz_output = ret

    rates_output_by_sat = {}
    rates_output_by_sat["obs"] = rates_output["obs"]
    rates_output_by_sat["xlink_update"] = rates_output["xlink_update"]
    rates_output_by_sat["xlink_rates_update"] = rates_output["xlink_rates_update"]
    rates_output_by_sat["xlnkrange_update"] = rates_output["xlnkrange_update"]
    rates_output_by_sat["gslink_update"] = rates_output["gslink_update"]
    rates_output_by_sat["gsaer_update"] = rates_output["gsaer_update"]
    rates_output_by_sat["gslink_rates_update"] = rates_output["gslink_rates_update"]

    viz_output_by_sat = {}
    viz_output_by_sat["t_o"] = viz_output["t_o"]
    viz_output_by_sat["o_locations"] = viz_output["o_locations"]
    viz_output_by_sat["t_d"] = viz_output["t_d"]
    viz_output_by_sat["d_part"] = viz_output["d_part"]
    viz_output_by_sat["t_x"] = viz_output["t_x"]
    viz_output_by_sat["x_part"] = viz_output["x_part"]
    viz_output_by_sat["dlnk_rate_history"] = viz_output["dlnk_rate_history"]
    viz_output_by_sat["xlnk_rate_history"] = viz_output["xlnk_rate_history"]

    return rates_output_by_sat, viz_output_by_sat
