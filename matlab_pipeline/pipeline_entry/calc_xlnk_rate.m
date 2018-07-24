function [xlnk_rate_Mbps] = calc_xlnk_rate(range_km,P_tx_xlnk,HPBW_xlnk,pointing_error,range_m,comm_type_xlnk_indx,params_struct,comm_type_xlnk_builtin,xlnk_table_interpolation,xlnk_table_ranges,xlnk_table_range_rates)

    % if using built-in link budget
    if comm_type_xlnk_builtin
        range_m = range_km*1000; % in meters
        %  note that this immediately gets converted to Mbps
        xlnk_rate_Mbps = crosslink_capacity(P_tx_xlnk,HPBW_xlnk,pointing_error,range_m,comm_type_xlnk_indx,params_struct)/1000000;
    % if using a provided lookup table
    else
        if strcmp (xlnk_table_interpolation,'floor') 
            %  get all ranges in the table less than our current range, then take the last of those ("flooring" the current range)
            rate_entry_indx_unfilt = find(range_km-xlnk_table_ranges>=0);
            rate_entry_indx = rate_entry_indx_unfilt(end);
        end

        xlnk_rate_Mbps = xlnk_table_range_rates(rate_entry_indx,2);

    end
    
end