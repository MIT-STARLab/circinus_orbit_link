def get_sat_dlnk_params(comm_type_dlnk, comm_type_dlnk_builtin,
                        comm_types_builtin, params):
    comm_type_dlnk_indx = -9
    dlnk_table_range_rates = []
    dlnk_table_ranges = []
    dlnk_table_interpolation = 'floor'

    # Comment ported from matlab:
    # determine the built-in index if specified, otherwise grab the rates lookup table

    if comm_type_dlnk_builtin:
        comm_type_dlnk_indx = comm_types_builtin[comm_type_dlnk]
    else:
        # Comment ported from matlab:
        # note that the .(comm_type_dlnk) below actually looks up the field name in the structure from a string

        # TODO: Pass the proper key in instead of this hack
        dlnk_dict = {}
        keys = {
            "range_rates_table",
            "range_units",
            "rates_units",
            "interpolation_method"
        }
        for val in params["lookup_params"]["dlnk_range_rates"].values():
            try:
                if keys.issubset(val.keys()):
                    dlnk_dict = val
            except AttributeError:
                pass  # This val was not a dict
        assert len(dlnk_dict) > 0

        dlnk_table_range_rates = dlnk_dict["range_rates_table"]

        dlnk_table_ranges = [i[0] for i in dlnk_table_range_rates]

        units_range = dlnk_dict["range_units"]
        units_rates = dlnk_dict["rates_units"]
        dlnk_table_interpolation = dlnk_dict["interpolation_method"]

        # TODO: use pint for this
        if units_range != 'km':
            raise ValueError(
                "{} range units are not implemented for xlnk range rates table".format(units_range))
        if units_rates != 'Mbps':
            raise ValueError(
                '{} rates units are not implemented for xlnk range rates table'.format(units_rates))
        if dlnk_table_interpolation != 'floor':
            raise ValueError('{} interpolation method is not implemented for xlnk range rates table'.format(
                xlnk_table_interpolation))

    return comm_type_dlnk_indx, dlnk_table_range_rates, dlnk_table_ranges, dlnk_table_interpolation
