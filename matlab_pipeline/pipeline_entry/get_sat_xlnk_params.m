function [comm_type_xlnk_indx,xlnk_table_range_rates,xlnk_table_ranges,xlnk_table_interpolation] = get_sat_xlnk_params(comm_type_xlnk,comm_type_xlnk_builtin,comm_types_builtin,params)

    comm_type_xlnk_indx = -9;
    xlnk_table_range_rates = [];
    xlnk_table_ranges = [];
    xlnk_table_interpolation = 'floor';

    if comm_type_xlnk_builtin
        comm_type_xlnk_indx = comm_types_builtin(comm_type_xlnk);  
    else
        % note that the .(comm_type_xlnk) below actually looks up the field name in the structure from a string
        xlnk_table_range_rates = params.lookup_params.xlnk_range_rates.(comm_type_xlnk).range_rates_table;
        xlnk_table_ranges = xlnk_table_range_rates(:,1);
        units_range = params.lookup_params.xlnk_range_rates.(comm_type_xlnk).range_units;
        units_rates = params.lookup_params.xlnk_range_rates.(comm_type_xlnk).rates_units;
        xlnk_table_interpolation = params.lookup_params.xlnk_range_rates.(comm_type_xlnk).interpolation_method;

        if ~strcmp(units_range,'km')
            error(disp('%s range units are not implemented for xlnk range rates table',units_range));
        end
        if ~strcmp(units_rates,'Mbps')
            error(disp('%s rates units are not implemented for xlnk range rates table',units_rates));
        end
        if ~strcmp(xlnk_table_interpolation,'floor')
            error(disp('%s interpolation method is not implemented for xlnk range rates table',xlnk_table_interpolation));
        end
    end
end