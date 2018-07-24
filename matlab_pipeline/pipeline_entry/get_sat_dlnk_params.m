function [comm_type_dlnk_indx,dlnk_table_range_rates,dlnk_table_ranges,dlnk_table_interpolation] = get_sat_dlnk_params(comm_type_dlnk,comm_type_dlnk_builtin,comm_types_builtin,params)
    comm_type_dlnk_indx = -9;
    dlnk_table_range_rates = [];
    dlnk_table_ranges = [];
    dlnk_table_interpolation = 'floor';

    %  determine the built-in index if specified, otherwise grab the rates lookup table
    if comm_type_dlnk_builtin
        comm_type_dlnk_indx = comm_types_builtin(comm_type_dlnk);  
    else
        % note that the .(comm_type_dlnk) below actually looks up the field name in the structure from a string
        dlnk_table_range_rates = params.lookup_params.dlnk_range_rates.(comm_type_dlnk).range_rates_table;
        dlnk_table_ranges = dlnk_table_range_rates(:,1);
        units_range = params.lookup_params.dlnk_range_rates.(comm_type_dlnk).range_units;
        units_rates = params.lookup_params.dlnk_range_rates.(comm_type_dlnk).rates_units;
        dlnk_table_interpolation = params.lookup_params.dlnk_range_rates.(comm_type_dlnk).interpolation_method;

        if ~strcmp(units_range,'km')
            error(disp('%s range units are not implemented for dlnk range rates table',units_range));
        end
        if ~strcmp(units_rates,'Mbps')
            error(disp('%s rates units are not implemented for dlnk range rates table',units_rates));
        end
        if ~strcmp(dlnk_table_interpolation,'floor')
            error(disp('%s interpolation method is not implemented for dlnk range rates table',dlnk_table_interpolation));
        end
    end
end