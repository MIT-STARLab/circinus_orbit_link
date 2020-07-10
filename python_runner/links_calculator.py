from datetime import datetime
from math import pi
from get_sat_xlnk_params import get_sat_xlnk_params
from get_sat_dlnk_params import get_sat_dlnk_params
import numpy as np

SEC_PER_DAY = 60 * 60 * 24
MICROSEC_PER_DAY = 1e6 * SEC_PER_DAY

def from_iso_format(dtstr):
    # Based on 3.7 datetime module
    # TODO: Support higher precision than hour

    year = int(dtstr[0:4])
    if dtstr[4] != '-':
        raise ValueError('Invalid date separator: %s' % dtstr[4])

    month = int(dtstr[5:7])

    if dtstr[7] != '-':
        raise ValueError('Invalid date separator: %s' % dtstr[7])

    day = int(dtstr[8:10])

    if dtstr[10] != 'T':
        raise ValueError('Invalid date separator: %s' % dtstr[10])

    hour = int(dtstr[11:13])

    return datetime(year, month, day, hour)


def mjuliandate(date_str):
    """
    Replacement for matlab mjuliandate function

    :param date_str: a string describing the date in this format:
      YYYY-mm-ddTHH:MM:SS.SSSSSSZ
    :return: modified julian date
    """

    # According to Wikipedia, Modified Julian Date starts at 0h Nov 17, 1858
    mjd_start = datetime(1858, 11, 17)

    # Strip Z at end
    assert date_str[-1] == 'Z'
    date_str = date_str[:-1]
    date_in = from_iso_format(date_str)

    # Returns a time_delta object
    mjd = date_in - mjd_start

    return mjd.days + mjd.seconds/SEC_PER_DAY + mjd.microseconds/MICROSEC_PER_DAY


def zeros(*args):
    """Like numpy.zeros but for lists."""
    if len(args) == 0:
        raise ValueError('"zeros" requires at least one dimension')
    elif len(args) == 1:
        return [0] * args[0]
    else:
        out = []
        for _ in range(args[0]):
            out.append(zeros(*args[1:]))
        return out


def lists(*args):
    if len(args) == 0:
        return []
    else:
        out = []
        for _ in range(args[0]):
            out.append(lists(*args[1:]))
        return out


def find(arr, func):
    ret = []
    for i, val in enumerate(arr):
        if func(val):
            ret.append(i)
    return ret

def calc_xlnk_rate(range_km,P_tx_xlnk, HPBW_xlnk,pointing_error,
        range_m, comm_type_xlnk_indx, params_struct, comm_type_xlnk_builtin,
        xlnk_table_interpolation, xlnk_table_ranges, xlnk_table_range_rates):
    # Not moved in a separate file like in matlab because of find func

    # if using built-in link budget
    if comm_type_xlnk_builtin:
        raise NotImplementedError
        # NOTE: This section uses a function "crosslink_capacity" that I can't
        # find anywhere in the SPRINT repo--I think "comm_type_xlnk_builtin" is
        # a deprecated option. But this is what the code for this section looks
        # like in python:
        # range_m = range_km * 1000 # in meters -- this also being passed as a
        # a parameter is a sign of likely deprecated code
        # # note that this immediately gets converted to Mbps
        # xlnk_rate_Mbps = crosslink_capacity(P_tx_xlnk, HPBW_xlnk, pointing_error,
        #     range_m, comm_type_xlnk_indx, params_struct) / 1000000:
    else: # if using a provided lookup table
        if xlnk_table_interpolation == 'floor':
            # get all ranges in the table less than our current range, then take the last of those ("flooring" the current range)
            rate_entry_indx_unfilt = find(xlnk_table_ranges, lambda x:
                range_km - x >= 0)
            rate_entry_indx = rate_entry_indx_unfilt[-1]

        xlnk_rate_Mbps = xlnk_table_range_rates[rate_entry_indx][1]

    return xlnk_rate_Mbps


def links_calculator(accesses_data, params, verbose):
    """
    Parameters: accesses_data, params, verbose
    Returns:
    """

    if verbose:
      print("running links calculator")

    # Set up parameters
    # TODO: Figure out what accessees_data is, how it is best used, etc
    obs  = accesses_data["obs"]
    gslink  = accesses_data["gslink"]
    gsaer  = accesses_data["gsaer"]
    xlink  = accesses_data["xlink"]
    xrange  = accesses_data["xrange"]

    # TODO: This will probably involve some function from time,
    # hence the import at the top
    start_time_dt = params["scenario_start_utc"]

    # More parameter set up
    # TODO: Similar as above.
    xlnk_max_len_secs = params["general_link_params"]["xlnk_max_len_s"][0]
    dlnk_max_len_secs = params["general_link_params"]["dlnk_max_len_s"][0]
    xlnk_storage_decimation = params["general_link_params"]["xlnk_decimation"]
    dlnk_storage_decimation = params["general_link_params"]["dlnk_decimation"]

    num_sats = params["num_sats"]
    num_gs = params["num_gs"]
    num_targets = params["num_targets"]

    DTR = pi/180

    # NOTE: This seemes like an enum...we'll see if it gets used that way
    # TODO: If so, does this UHF-Cadet can be made to be 5 instead?
    comm_types_builtin = {}
    comm_types_builtin["Optical-general"] = 0
    comm_types_builtin["Ka-band-general"] = 1
    comm_types_builtin["X-band-general"] = 2
    comm_types_builtin["Optical-Sinclair"] = 3
    comm_types_builtin["Optical-OCSD"] = 4
    comm_types_builtin["UHF-Cadet"] = 6

    t_0 = mjuliandate(start_time_dt)

    # Output data
    rates_output = {}

    # Temporary containers
    gslink_rates = lists(num_sats, num_gs)
    xlink_rates = lists(num_sats, num_sats)
    xlink_directions = lists(num_sats, num_sats)

    # Data rate outputs
    xlink_update = lists(num_sats, num_sats)
    xlink_rates_update  = lists(num_sats, num_sats)
    xlnkrange_update = lists(num_sats, num_sats)

    gslink_update = lists(num_sats, num_gs)
    gsaer_update = lists(num_sats, num_gs)
    gslink_rates_update = lists(num_sats, num_gs)

    # visualization output
    t_o = [[] for _ in range(num_sats)]
    o_locations = [[] for _ in range(num_sats)]

    t_d = lists(num_sats)
    d_part = lists(num_sats)
    t_x = lists(num_sats)
    x_part = lists(num_sats)
    dlnk_rate_history = lists(num_sats, num_gs)
    xlnk_rate_history = lists(num_sats, num_sats)

    # Original comment: these structures indicate which downlink and cross-link pairs are to be discluded from processing in the code below. if a particular satellite ID is present as a field within the structure, then the list of disableded ground station or satellite IDs are strings within a cell array for that field.
    # My modified comment: If a satellite ID is in the keys of this dict, then the disabled ground station or satellite IDs are in a set that is the value
    # TODO: Ensure that the params follow this convention
    link_direc_disables_dlnk = params["link_disables"]["dlnk_direc_disabled_gs_indx_by_sat_indx"]
    link_direc_disables_xlnk = params["link_disables"]["xlnk_direc_disabled_xsat_indx_by_sat_indx"]

    gs_params = params["gs_params"]

    # Run
    for sat_num in range(num_sats):
        if verbose:
            print("sat num:", sat_num)

        # Check for disabled IDs
        link_direc_disables_dlnk_sat = set()
        link_direc_disables_xlnk_sat = set()

        try:
            link_direc_disables_dlnk_sat = link_direc_disables_dlnk[sat_num]
        except KeyError:
            pass

        try:
            link_direc_disables_xlnk_sat = link_direc_disables_xlnk[sat_num]
        except KeyError:
            pass

        # Get satellite parameters
        # These two are not actually used in currently run
        P_tx_dlnk_sat = params["sat_link_params"][sat_num]["dlnk_params"]["P_tx_W"]
        HPBW_dlnk_sat = params["sat_link_params"][sat_num]["dlnk_params"]["HPBW_rad"] / 1000

        P_tx_xlnk_sat = params["sat_link_params"][sat_num]["xlnk_params"]["P_tx_W"]
        HPBW_xlnk_sat = params["sat_link_params"][sat_num]["xlnk_params"]["HPBW_rad"] / 1000
        pointing_error_sat = params["sat_link_params"][sat_num]["pointing_error_deg"] * DTR

        # the comm type for dlnk, specifies which params to use
        comm_type_dlnk = params["sat_link_params"][sat_num]["dlnk_params"]["comm_type"]["name"]
        # specifies whether or not we're using a built-in model for the downlink link budget or a lookup table from the input params
        comm_type_dlnk_builtin = params["sat_link_params"][sat_num]["dlnk_params"]["comm_type"]["built_in"]
        comm_type_xlnk_sat = params["sat_link_params"][sat_num]["xlnk_params"]["comm_type"]["name"]
        comm_type_xlnk_sat_builtin = params["sat_link_params"][sat_num]["xlnk_params"]["comm_type"]["built_in"]

        [comm_type_dlnk_indx, dlnk_table_range_rates, dlnk_table_ranges,
            dlnk_table_interpolation] = get_sat_dlnk_params(
                comm_type_dlnk, comm_type_dlnk_builtin, comm_types_builtin, params)


        [comm_type_xlnk_sat_indx, xlnk_sat_table_range_rates, xlnk_sat_table_ranges,
            xlnk_sat_table_interpolation] = get_sat_xlnk_params(
                comm_type_xlnk_sat, comm_type_xlnk_sat_builtin, comm_types_builtin, params)

        # Handle observations - note this is only for viz in cesiumJS
        if verbose:
            print("obs")

        for targ_num in range(num_targets):
            num_obs = len(obs[sat_num][targ_num])

            for obs_num in range(num_obs):
                tstart = obs[sat_num][targ_num][obs_num][0]
                tstop = obs[sat_num][targ_num][obs_num][1]

                t_o[sat_num].append([tstart, tstop])
                o_locations[sat_num].append(int(targ_num))

        # Matlab comment:
        # have to sort t_o, because obs show times are not created uniquely for
        # each sat-obs target combination, instead they're all lumpted into a
        # single satellite (in VizInputsGenerator.py for cesiumJS input file)
        for sat_t_o in t_o:
            # Sort by start time (first element), push empty lists to the back
            sat_t_o.sort(key=lambda obs: obs[0])

        # Handle downlinks
        if verbose:
            print("dlnks")

        # Calculate data rates at each timestep of downlink and store in gslink_rates
        for gs_num in range(num_gs):
            comm_type_gs = gs_params['stations'][gs_num]["comm_type"]
            # check if the dlnk types are compatible (their names are equal).
            # if not, skip this gs
            # TODO: use enums instead
            if comm_type_dlnk != comm_type_gs:
                print("skipping gs index {} comm_type does not match".format(gs_num))
                continue

            # Check if this gs is disabled for the satellite
            if gs_num in link_direc_disables_dlnk_sat:
                continue

            gslink_rates[sat_num][gs_num] = zeros(len(gsaer[sat_num][gs_num]))
            num_dlnks = len(gslink[sat_num][gs_num])  # TODO: double check

            # COPIED TODO: add compatibility check between sat and gs comm types

            for dlnk_num in range(num_dlnks):
                num_timepoints = len(gsaer[sat_num][gs_num][dlnk_num])
                gslink_rates[sat_num][gs_num][dlnk_num] = zeros(num_timepoints, 2)

                # TODO: Probably wanna change i to t...
                for i in range(num_timepoints):
                    gslink_rates[sat_num][gs_num][dlnk_num][i][0] = gsaer[sat_num][gs_num][dlnk_num][i][0]  # store time
                    # elev = gsaer[sat_num][gs_num][dlnk_num][i][2]
                    range_km = gsaer[sat_num][gs_num][dlnk_num][i][3]

                    # if using built-in link budget
                    if comm_type_dlnk_builtin:
                        # NOTE: I think this is a deprecated code branch. I don't think
                        # comm_type_dlnk_builtin is ever true.
                        raise NotImplementedError

                        # What the approximate translated version would have looked like:
                        # (NOTE: downlink_capacity doesn't appear to exist anywhere, even
                        # when I grepped the whole SPRINT repo)
                        #
                        # range_m = range_km * 1000  # Convert to meters
                        # # note that this immediately gets converted to Mbps
                        # gslink_rates[sat_num][gs_num][dlnk_num][i][1] = downlink_capacity(P_tx_dlnk_sat, HPBW_dlnk_sat, pointing_error_sat, range_m, comm_type_dlnk_indx, elev, gs_num, params_struct)/1000000  # TODO: make this a real function
                    else:
                        if dlnk_table_interpolation == "floor":
                            # get all ranges in the table less than our current range, then take the last of those ("flooring" the current range)
                            rate_entry_indx_unfilt = find(dlnk_table_ranges,
                                lambda x: range_km - x >= 0)
                            # NOTE: This assumes they start in ascending order
                            # This is guaranteed with the current find implementation
                            rate_entry_indx = rate_entry_indx_unfilt[-1]

                        gslink_rates[sat_num][gs_num][dlnk_num][i][1] = dlnk_table_range_rates[rate_entry_indx][1]

                    assert np.imag(gslink_rates[sat_num][gs_num][dlnk_num][i][1]) == 0

        # Copy of old TODO:
        # TODO: probably ought to trim off any zero downlink data rate time
        # from gslink, gsaer, gslink_rates...left a warning in the code below
        # to catch and warn about this situation...

        # divide downlink windows into smaller time periods to make them more
        # manageable, as well as decimate the gsaer and data rates storage
        dlnk_max_len_dayf = dlnk_max_len_secs / 86400

        for gs_num in range(num_gs):
            comm_type_gs = gs_params['stations'][gs_num]["comm_type"]
            # check if the dlnk types are compatible (their names are equal).  if not, skip this gs
            # TODO: use enums instead
            if comm_type_dlnk != comm_type_gs:
                continue

            # Check if this gs is disabled for the satellite
            if gs_num in link_direc_disables_dlnk_sat:
                continue

            gslink_update[sat_num][gs_num] = []

            num_dlnks = len(gslink[sat_num][gs_num])

            gsaer_update[sat_num][gs_num] = []
            gslink_rates_update[sat_num][gs_num] = []
            for dlnk_num in range(num_dlnks):
                if gslink_rates[sat_num][gs_num][dlnk_num] == 0:
                    # TODO: actually trim this
                    # NOTE: matlab code also had a warning for this at a different
                    # line
                    print("Skipping over empty dlnks_rate at dlnk_num={}".format(
                        dlnk_num))
                    continue

                dlnk_len =  gslink[sat_num][gs_num][dlnk_num][2]

                dlnk_aer_indx = 0
                dlnk_aer_indx_start = 0

                # NOTE: With this t, the following loop is unreachable, but
                # it's fine (see comment in matlab code)
                t = dlnk_max_len_dayf

                dlnk_start = gslink[sat_num][gs_num][dlnk_num][0]
                dlnk_end = gslink[sat_num][gs_num][dlnk_num][1]

                last_t = 0

                # loop through as long as end time of slice is less than
                # the end of the full xlnk access
                # while dlnk_start + t < dlnk_end:
                #     print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                #     # NOTE: In the test inputs from Bobby, thie code was never reached
                #     # (t was too short)

                #     # TODO: Fix this shape
                #     decimation = dlnk_storage_decimation[0]

                #     # figure out up to what index in xlnk_xrange (same for
                #     # xlink_rates_update[sat_num][xsat_num][xlnk_num]) we need to slice
                #     # TODO: There is probably a better built in function for this
                #     new_end = dlnk_start + t
                #     dlnk_aer_indx_start = dlnk_aer_indx
                #     while gsaer[sat_num][gs_num][dlnk_num][dlnk_aer_indx][0] <= new_end:
                #         dlnk_aer_indx += 1
                #         if dlnk_aer_indx == len(gsaer[sat_num][gs_num][dlnk_num]) - 1:
                #             break

                #     # before we add the new slice check if there are enough indices
                #     # for it to be relevant to add
                #     if (dlnk_aer_indx - dlnk_aer_indx_start <
                #             decimation):
                #         break

                #     gslink_update[sat_num][gs_num].append([dlnk_start + last_t, dlnk_start + t])

                #     # Decimate
                #     dlnk_aer_slice = gsaer[sat_num][gs_num][dlnk_num][dlnk_aer_indx_start:dlnk_aer_indx:
                #         decimation]
                #     dlnk_rates_slice = gsaer_update[sat_num][gs_num][dlnk_aer_indx_start:dlnk_aer_indx:
                #         decimation]

                #     gsaer_update[sat_num][gs_num].append(dlnk_aer_slice)
                #     gslink_rates_update[sat_num][gs_num].append(dlnk_rates_slice)

                #     last_t = t
                #     t = t + dlnk_max_len_dayf
                # handle last slice
                gslink_update[sat_num][gs_num].append([dlnk_start + last_t, dlnk_end])
                dlnk_aer_indx_start = dlnk_aer_indx
                dlnk_aer_slice = gsaer[sat_num][gs_num][dlnk_num][dlnk_aer_indx_start::dlnk_storage_decimation[0]]
                gsaer_update[sat_num][gs_num].append(dlnk_aer_slice)
                gslink_rates_update[sat_num][gs_num].append(
                    gslink_rates[sat_num][gs_num][dlnk_num][dlnk_aer_indx_start::dlnk_storage_decimation[0]])

        # Comment from matlab:
        # format dlnk history for input to CesiumJS. Note of course that this
        # is super inefficient, but this code only runs once, and splitting all
        # this code out helps preserve my sanity
        for gs_num in range(num_gs):
            comm_type_gs =  gs_params['stations'][gs_num]['comm_type']

            if comm_type_dlnk != comm_type_gs:
                continue

            # Skip disabled gs for this sattellite
            if any(gs_num == i for i in link_direc_disables_dlnk_sat):
                continue

            num_dlnks = len(gslink_update[sat_num][gs_num])

            for dlnk_num in range(num_dlnks):
                num_timepoints = len(gsaer_update[sat_num][gs_num][dlnk_num])

                # Matlab comment:
                # store in dlnk history for visualizing in cesiumJS

                # Matlab comment:
                # Need to use dlnk for this point because times in gsaer are
                # decimated
                tstart = gslink_update[sat_num][gs_num][dlnk_num][0]
                # NOTE: This was originally index (1, 1) which is allowed in matlab
                # even for 1D arrays

                # 0 rate right before start of downlink
                new_time = (tstart - t_0) * 86400 - 0.5
                dlnk_rate_history[sat_num][gs_num].append([new_time - 0.1, 0])
                # (minus 0.1 for ~impulsive transition)
                last_history_val = 0

                for i in range(num_timepoints):
                    # Matlab comment:
                    # decimate dlnk rate storage history, to maintain some
                    # sanity
                    # note decimating on top of decimation already done above

                    # NOTE: The line before had a comment inserted before ==,
                    # meaning there was either a bug here or some code which was
                    # never reached
                    if i % 2 == 1:
                        new_time = (gsaer_update[sat_num][gs_num][dlnk_num][i][0] - t_0) * 86400

                        # add point for rate we're transitioning from (so transition is
                        # impulsive)
                        dlnk_rate_history[sat_num][gs_num].append([new_time, last_history_val])

                        dlnk_rate_history[sat_num][gs_num].append([
                            new_time + 0.1,
                            gslink_rates_update[sat_num][gs_num][dlnk_num][i][1]
                        ])

                        # add point for date rate during pass (plus 0.1 for ~impulsive
                        # transition)
                        last_history_val = gslink_rates_update[sat_num][gs_num][dlnk_num][i][1]

                # Matlab comment:
                # Need to use dlnk for this point because times in gsaer are
                # decimated
                tstop = gslink_update[sat_num][gs_num][dlnk_num][1]
                new_time = (tstop - t_0) * 86400

                dlnk_rate_history[sat_num][gs_num].append([new_time + 0.2, last_history_val])
                # Matlab comment:
                # make sure to add in a point at end of pass, so final rate dropoff
                # is "impulsive". 0.2 so it'll come after any point added above
                # (could possibly conflict with +0.1 above)
                dlnk_rate_history[sat_num][gs_num].append([new_time + 0.3, 0])
                # now zero it (plus 0.3 so it'll be, again, after any points above)

                t_d[sat_num].append([tstart, tstop])
                # explicitly save partner variables as integers
                d_part[sat_num].append(int(gs_num))  # NOTE: was uint16
        # Matlab comment:
        # note: don't need to sort t_d, cesium generation code can
        # handle it, because downlink show times are created separately for
        # each sat-gs combination in VizInputsGenerator.py for cesiumJS input
        # file


        # Matlab comment:
        # Handle crosslinks

        if verbose:
            print('xlnks')

        skip_xlnk = zeros(num_sats, num_sats)
        # Matlab comment:
        # calculate datarates at each timestep during crosslinks - store in new
        # cell array xlink_rates
        print("Sat num: {}".format(sat_num+1, num_sats))
        # TODO: should probably remove this print at the end
        for xsat_num in range(sat_num+1, num_sats):
            # Start from sat_num+1 because xlink is symmetric

            comm_type_xlnk_xsat = params['sat_link_params'][1]['xlnk_params'][
                'comm_type']['name']

            # TODO; What is the point of these?
            skip_xlnk_sat_xsat = False
            skip_xlnk_xsat_sat = False

            # Matlab comment:
            # check if the crosslink types are compatible (their names are equal).
            # if not, skip this satellite
            # TODO: If different technologies are desired for TX and RX, need to
            # update this check here. for now are assuming that both satellites
            # must have the same technology
            if comm_type_xlnk_sat != comm_type_xlnk_xsat:
                if verbose:
                    print('skipping xsat python indx', xsat_num)
                continue

            # Matlab comment:
            # determine if the cross-link is explicitly disabled in either
            # direction. if yes, then we'll be skipping a direction
            link_direc_disables_xlnk_xsat = set()
            try:
                link_direc_disables_xlnk_xsat = link_direc_disables_xlnk[xsat_num]
            except KeyError:
                pass

            # check if the xsat index is within the disabled list for the satellite. If yes then skip this satellite
            skip_xlnk_sat_xsat = skip_xlnk_sat_xsat or xsat_num in \
                link_direc_disables_xlnk_sat

            # check if the xsat index is within the disabled list for the
            # satellite. If yes then skip this satellite
            skip_xlnk_xsat_sat = skip_xlnk_xsat_sat or sat_num in \
                link_direc_disables_xlnk_xsat

            if skip_xlnk_sat_xsat:
                if verbose:
                    print('skipping xlnk', sat_num, '->', xsat_num)

            if skip_xlnk_xsat_sat:
                if verbose:
                    print('skipping xlnk', xsat_num, '->', sat_num)

            # Matlab comment:
            # only don't procede through the rest of this iteration if we are
            # skipping both directions in the cross-link
            if skip_xlnk_sat_xsat and skip_xlnk_xsat_sat:
                skip_xlnk[sat_num][xsat_num] = 1
                continue


            P_tx_xlnk_xsat = params['sat_link_params'][xsat_num]['xlnk_params']['P_tx_W']
            HPBW_xlnk_xsat = params['sat_link_params'][xsat_num]['xlnk_params']['HPBW_rad']/1000
            pointing_error_xsat = params['sat_link_params'][xsat_num]['pointing_error_deg'] * DTR

            comm_type_xlnk_xsat_builtin = params['sat_link_params'][xsat_num]['xlnk_params'][
                'comm_type']['built_in']
            [comm_type_xlnk_xsat_indx, xlnk_xsat_table_range_rates,
                xlnk_xsat_table_ranges, xlnk_xsat_table_interpolation] = \
                    get_sat_xlnk_params(comm_type_xlnk_xsat,
                    comm_type_xlnk_xsat_builtin, comm_types_builtin,params)

            # Matlab comment:
            # record if they are the same params, in which case the calculated data
            # rates will be the same so we can save some processing
            same_params_both_dirs = (P_tx_xlnk_sat == P_tx_xlnk_xsat) and \
                (HPBW_xlnk_sat == HPBW_xlnk_xsat) and (pointing_error_sat == pointing_error_xsat) and \
                (comm_type_xlnk_sat == comm_type_xlnk_xsat)

            xlink_rates[sat_num][xsat_num] = zeros(len(xrange[sat_num][xsat_num]))
            xlink_directions[sat_num][xsat_num] = zeros(len(xrange[sat_num][xsat_num]))

            num_xlnks = len(xlink[sat_num][xsat_num])

            print('Num xlnks:', num_xlnks)

            for xlnk_num in range(num_xlnks):
                num_timepoints = len(xrange[sat_num][xsat_num][xlnk_num])
                # for each row, first index is timestamp, second index is data rate
                # from satellite to cross-link satellite, and third index is data
                # rate from cross-link satellite to satellite
                xlnk_rates = zeros(num_timepoints, 3)

                for i in range(num_timepoints):
                    xlnk_rates[i][0] = xrange[sat_num][xsat_num][xlnk_num][i][0]
                    range_km = xrange[sat_num][xsat_num][xlnk_num][i][1]

                    if not comm_type_xlnk_sat_builtin and \
                            xlnk_sat_table_interpolation == 'floor':
                        # get all ranges in the table less than our current range,
                        # then take the last of those ("flooring" the current
                        # range)
                        rate_entry_indx_unfilt = find(xlnk_sat_table_ranges,
                            lambda x: range_km - x >= 0)
                        rate_entry_indx_sat = rate_entry_indx_unfilt[-1]

                    if skip_xlnk_sat_xsat:
                        xlnk_rates[i][1] = 0
                    else:
                        if comm_type_xlnk_sat_builtin:
                            # I think the above bool is never true ^
                            # 1) range_m is not defined in the matlab code, I had
                            #    to go back and do the conversion here
                            # 2) comm_type_xlnk_sat_builtin also seems like a
                            #    not used case inside calc_xlnk_rate
                            # 3) Also I just realized that in the matlab code,
                            #    at the top it says that builtin comm types are
                            #    not being used
                            # 4) There are other undefined parameters here

                            raise NotImplementedError

                            # This is what the python function call would look
                            # like at least
                            # range_m = range_km * 1000
                            # xlnk_rates[i][1] = calc_xlnk_rate(
                            #     range_km,P_tx_xlnk_sat,
                            #     HPBW_xlnk_sat,
                            #     pointing_error_sat,
                            #     range_m,
                            #     comm_type_xlnk_sat_indx,
                            #     params_struct,
                            #     comm_type_xlnk_sat_builtin[0],
                            #     xlnk_sat_table_interpolation,
                            #     xlnk_sat_table_ranges,
                            #     xlnk_sat_table_range_rates)
                        else:
                            xlnk_rates[i][1] = xlnk_sat_table_range_rates[
                                rate_entry_indx_sat][1]

                    # Matlab comment:
                    # FML, this is so ugly...
                    if skip_xlnk_xsat_sat:
                        xlnk_rates[i][2] = 0
                    else:
                        if same_params_both_dirs and not skip_xlnk_sat_xsat:
                            xlnk_rates[i][2] = xlnk_rates[i][1]
                        else:
                            if comm_type_xlnk_xsat_builtin:
                                # See previous comment for comm_type_xlnk_xsat_builtin
                                raise NotImplementedError
                                # xlnk_rates[i][2] = calc_xlnk_rate(range_km,P_tx_xlnk_xsat,HPBW_xlnk_xsat,pointing_error_xsat,range_m,comm_type_xlnk_xsat_indx,params_struct,comm_type_xlnk_xsat_builtin,xlnk_xsat_table_interpolation,xlnk_xsat_table_ranges,xlnk_xsat_table_range_rates)
                            else:
                                xlnk_rates[i][2] = xlnk_xsat_table_range_rates[
                                    rate_entry_indx_xsat][1]

                    # Matlab comment:
                    # leaving in this, emily's old diagnostic code, just in
                    # case
                    if np.imag(xlnk_rates[i][1]) > 0:
                        print('imaginary part in xlnk rate...setting to 73 mbps for now')

                        xlnk_rates[i][1] = 73.1133803921651
                    if np.imag(xlnk_rates[i][2]) > 0:
                        print('imaginary part in xlnk rate...setting to 73 mbps for now')

                        xlnk_rates[i][2] = 73.1133803921651
                xlink_directions[sat_num][xsat_num][xlnk_num] = zeros(3)
                if not skip_xlnk_sat_xsat:
                    # record a "true" for the first satellite transmitting
                    xlink_directions[sat_num][xsat_num][xlnk_num][0] = 1
                if not skip_xlnk_xsat_sat:
                    # record a "true" for the second satellite transmitting
                    xlink_directions[sat_num][xsat_num][xlnk_num][1] = 1
                if same_params_both_dirs:
                    # record a "true" for the satellites transmitting at the same rate
                    xlink_directions[sat_num][xsat_num][xlnk_num][2] = 1

                xlink_rates[sat_num][xsat_num][xlnk_num] = xlnk_rates
        # Matlab comment:
        # divide crosslink windows into smaller time periods to make them more
        # manageable, as well as decimate the xrange and data rates storage
        # TODO: this seems like overlapping code with downlink
        xlnk_max_len_dayf = xlnk_max_len_secs / 86400
        for xsat_num in range(sat_num+1, num_sats):
            # want to start from sat_num+1 because xlink is symmetric:

            if skip_xlnk[sat_num][xsat_num]:
                continue

            # TODO: xrange shadows a builtin python funbc name, this should be changed
            num_xlnks = len(xlink[sat_num][xsat_num])

            xlink_update[sat_num][xsat_num] = []
            xlnks_xrange_update = []
            xlink_rates_update[sat_num][xsat_num] = []

            for xlnk_num in range(num_xlnks):
                xlnk_len = xlink[sat_num][xsat_num][xlnk_num][2]

                # Matlab comment:
                # if xlnk_len > xlnk_max_len_secs  # comment out because we:
                # want to run this code regardless - it just seemed like the
                # easiest way to include decimation of xrange, rates
                # grab these for slicing purposes
                xlnk_rates = xlink_rates[sat_num][xsat_num][xlnk_num]
                xlnk_xrange_indx = 0
                xlnk_xrange_indx_start = 0

                t = xlnk_max_len_dayf

                xlnk_start = xlink[sat_num][xsat_num][xlnk_num][0]
                xlnk_end = xlink[sat_num][xsat_num][xlnk_num][1]

                last_t = 0

                # Matlab comment:
                # loop through as long as end time of slice is less than
                # the end of the full xlnk access
                while xlnk_start + t < xlnk_end:
                    # Matlab comment:
                    # figure out up to what index in xrange[sat_num][xsat_num][xlnk_num] (same for
                    # xlnk_rates) we need to slice
                    new_end = xlnk_start + t
                    xlnk_xrange_indx_start = xlnk_xrange_indx
                    while xrange[sat_num][xsat_num][xlnk_num][xlnk_xrange_indx][0] <= new_end:
                        xlnk_xrange_indx += 1

                    decimation = xlnk_storage_decimation[0]

                    # Matlab comment:
                    # before we add the new slice, check if there are enough indices for it to be relevant to add]
                    if (xlnk_xrange_indx - xlnk_xrange_indx_start < decimation):
                        break

                    # Matlab comment:
                    # new sliced access window
                    # combine with the cross-linked directions array
                    xlink_update[sat_num][xsat_num].append([xlnk_start + last_t, xlnk_start + t] +
                        xlink_directions[sat_num][xsat_num][xlnk_num])

                    # Matlab comment:
                    # slice it. Note the -1 so we don't go past the end of
                    # the slice
                    # note the decimation factor used, decimation
                    xlnk_xrange_slice = xrange[sat_num][xsat_num][xlnk_num][xlnk_xrange_indx_start:
                        xlnk_xrange_indx:decimation]
                    xlnk_rates_slice = xlink_rates[sat_num][xsat_num][xlnk_num][xlnk_xrange_indx_start:
                        xlnk_xrange_indx:decimation]

                    xlnks_xrange_update.append(xlnk_xrange_slice)
                    xlink_rates_update[sat_num][xsat_num].append(xlnk_rates_slice)

                    last_t = t
                    t = t + xlnk_max_len_dayf

                # Matlab comment:
                # handle last slice
                xlink_update[sat_num][xsat_num].append([xlnk_start + last_t, xlnk_end] + xlink_directions[sat_num][xsat_num][xlnk_num])
                xlnk_xrange_indx_start = xlnk_xrange_indx
                xlnk_xrange_slice = xrange[sat_num][xsat_num][xlnk_num][
                    xlnk_xrange_indx_start::decimation]
                xlnk_rates_slice = xlink_rates[sat_num][xsat_num][xlnk_num][
                    xlnk_xrange_indx_start::decimation]
                xlnks_xrange_update.append(xlnk_xrange_slice)
                xlink_rates_update[sat_num][xsat_num].append(xlnk_rates_slice)

                # Matlab comment:
                # else # if we don't need to slice because xlnk isn't too long
                #     xlink_update[sat_num][xsat_num] = [xlink_update[sat_num][xsat_num]; xlnk]
                #     xlnks_xrange_update = [xlnks_xrange_update; xrange[sat_num][xsat_num][xlnk_num]]
                #     xlink_rates_update[sat_num][xsat_num] = [xlink_rates_update[sat_num][xsat_num]; xlink_rates[sat_num][xsat_num][xlnk_num]]
                # end

            xlnkrange_update[sat_num][xsat_num] = xlnks_xrange_update

        # Matlab comment:
        # TODO: this code currently ignores the fact the update to cross-link
        # rates that allows the two directions between satellitesto be
        # nonsymmetric. should update this to include it

        # Matlab comment:
        # format xlnk history for input to CesiumJS. Note of course that this
        # is super inefficient, but this code only runs once, and splitting all
        # this code out helps preserve my sanity
        for xsat_num in range(sat_num+1, num_sats):
            # want to start from sat_num+1 because xlink is symmetric:
            if skip_xlnk[sat_num][xsat_num]:
                continue

            num_xlnks = len(xlink_update[sat_num][xsat_num])

            for xlnk_num in range(num_xlnks):
                num_timepoints = len(xlnkrange_update[sat_num][xsat_num][xlnk_num])

                # Matlab comment:
                # store in xlnk history for visualizing in cesiumJS
                tstart = xlink_update[sat_num][xsat_num][xlnk_num][0]
                # Need to use xlnk for this point because times in xlnk_xrange are decimated

                new_time = (tstart - t_0) * 86400 - 0.5
                # 0 rate right before start of crosslink

                xlnk_rate_history[sat_num][xsat_num].append([new_time - 0.1, 0])
                # (minus 0.1 for ~impulsive transition)
                last_history_val = 0

                for i in range(num_timepoints):
                    if verbose:
                        if xlink_rates_update[sat_num][xsat_num][xlnk_num][i][1] == 0:
                            print('warning, xlink_rates_update[sat_num][xsat_num][xlnk_num][i][1] is zero. Perhaps this point should be trimmed...')

                    # Matlab comment:
                    # decimate xlnk rate storage history, to maintain some
                    # sanity
                    # note decimating on top of decimation already done above
                    if (i % 2) == 0:
                        new_time = (xlnkrange_update[sat_num][xsat_num][xlnk_num][i][0] - t_0) * 86400
                        # use time from range data struct because that captures
                        # every point in the window
                        xlnk_rate_history[sat_num][xsat_num].append([new_time, last_history_val])
                        # add point for rate we're transitioning from
                        # (so transition is impulsive)
                        xlnk_rate_history[sat_num][xsat_num].append(
                                [new_time + 0.1, xlink_rates_update[sat_num][xsat_num][xlnk_num][i][1]])
                            # add point for date rate during pass
                            # (plus 0.1 for ~impulsive transition)
                        last_history_val = xlink_rates_update[sat_num][xsat_num][xlnk_num][i][1]
                tstop = xlink_update[sat_num][xsat_num][xlnk_num][1]
                # Need to use xlnk for this point because times in
                # xlnk_xrange are decimated
                new_time = (tstop - t_0) * 86400
                xlnk_rate_history[sat_num][xsat_num].append(
                        [new_time + 0.2, last_history_val])
                    # make sure to add in a point at end of pass, so final rate
                    # dropoff is "impulsive". 0.2 so it'll come after any point
                    # added above (could possibly conflict with +0.1 above)
                xlnk_rate_history[sat_num][xsat_num].append(
                        [new_time + 0.3, 0])
                    # now zero it (plus 0.3 so it'll be, again, after any points above)


                t_x[sat_num].append([tstart, tstop])
                # explicitly save partner variables as integers
                x_part[sat_num].append(int(xsat_num)) # was uint16

        # Matlab comment:
        # note: don't need to sort t_x, cesium generation code can
        # handle it, because downlink show times are created separately for
        # each sat-gs combination in VizInputsGenerator.py for cesiumJS input
        # file


    # Matlab comment:
    ## Save output
    rates_output["obs"] = obs
    rates_output["xlink_update"] = xlink_update
    rates_output["xlink_rates_update"] = xlink_rates_update
    rates_output["xlnkrange_update"] = xlnkrange_update
    rates_output["gslink_update"] = gslink_update
    rates_output["gsaer_update"] = gsaer_update
    rates_output["gslink_rates_update"] = gslink_rates_update


    # Matlab comment:
    # save outputs for cesiumJS viz
    viz_output = {}
    viz_output["t_o"] = t_o
    viz_output["o_locations"] = o_locations
    viz_output["t_d"] = t_d
    viz_output["d_part"] = d_part
    viz_output["t_x"] = t_x
    viz_output["x_part"] = x_part
    viz_output["dlnk_rate_history"] = dlnk_rate_history
    viz_output["xlnk_rate_history"] = xlnk_rate_history


    # Matlab comment:
    # q_o_sizes_history
    # gs_availability_windows
    # batt_stored_history
    # t_eclipse

    return rates_output, viz_output
