def get_sat_xlnk_params(comm_type_xlnk, comm_type_xlnk_builtin,
                        comm_types_builtin, params):
    comm_type_xlnk_indx = -9
    xlnk_table_range_rates = []
    xlnk_table_ranges = []
    xlnk_table_interpolation = 'floor'

    # Comment ported from matlab:
    # determine the built-in index if specified, otherwise grab the rates lookup table

    if comm_type_xlnk_builtin:
        comm_type_xlnk_indx = comm_types_builtin["comm_type_xlnk"]
    else:
        # Comment ported from matlab:
        # note that the .(comm_type_xlnk) below actually looks up the field name in the structure from a string

        # TODO: Pass the proper key in instead of this hack
        xlnk_dict = {}
        keys = {
            "range_rates_table",
            "range_units",
            "rates_units",
            "interpolation_method"
        }
        for val in params["lookup_params"]["xlnk_range_rates"].values():
            try:
                if keys.issubset(val.keys()):
                    xlnk_dict = val
            except AttributeError:
                pass  # This val was not a dict
        assert len(xlnk_dict) > 0

        xlnk_table_range_rates = xlnk_dict["range_rates_table"]

        xlnk_table_ranges = [i[0] for i in xlnk_table_range_rates]

        units_range = xlnk_dict["range_units"]
        units_rates = xlnk_dict["rates_units"]
        xlnk_table_interpolation = xlnk_dict["interpolation_method"]

        # TODO: use pint for this
        if units_range != 'km':
            raise ValueError(
                "{} range units are not implemented for xlnk range rates table".format(units_range))
        if units_rates != 'Mbps':
            raise ValueError(
                '{} rates units are not implemented for xlnk range rates table'.format(units_rates))
        if xlnk_table_interpolation != 'floor':
            raise ValueError('{} interpolation method is not implemented for xlnk range rates table'.format(
                xlnk_table_interpolation))

    return comm_type_xlnk_indx, xlnk_table_range_rates, xlnk_table_ranges, xlnk_table_interpolation
