% Author: Kit Kennedy    
%  data rate calculations, as well as outputs for visualization engine

function [rates_output_by_sat,viz_output_by_sat] = links_wrapper(...
        accesses_data_by_sat,...
        params)

if params.verbose,
    disp('running links wrapper');
end

% find base directory for matlab pipeline code
file_dir = fileparts(mfilename('fullpath'));
base_directory = strcat(file_dir,'/..');  % matlab base

% add path to other code we'll be using
addpath(strcat(base_directory,'/code'));
addpath(strcat(base_directory,'/matlab_tools'));

%  reformat inputs
accesses_data = struct;
accesses_data.obs = reshape_inputs(accesses_data_by_sat.obs,params.num_targets);
accesses_data.gslink = reshape_inputs(accesses_data_by_sat.gslink,params.num_gs);
accesses_data.gsaer = reshape_inputs(accesses_data_by_sat.gsaer,params.num_gs);
accesses_data.xlink = reshape_inputs(accesses_data_by_sat.xlink,params.num_sats);
accesses_data.xrange = reshape_inputs(accesses_data_by_sat.xrange,params.num_sats);

% save('hey.mat','accesses_data','params','base_directory','params')

%  run the links calculator
[rates_output,viz_output] = links_calculator(...
        accesses_data,...
        params,...
        base_directory,...
        params.verbose);

%  reformat outputs
trim_empty = false;
rates_output_by_sat.obs = reshape_outputs(rates_output.obs,trim_empty);
rates_output_by_sat.xlink_update = reshape_outputs(rates_output.xlink_update,trim_empty);
rates_output_by_sat.xlink_rates_update = reshape_outputs(rates_output.xlink_rates_update,trim_empty);
rates_output_by_sat.xlnkrange_update = reshape_outputs(rates_output.xlnkrange_update,trim_empty);
rates_output_by_sat.gslink_update = reshape_outputs(rates_output.gslink_update,trim_empty);
rates_output_by_sat.gsaer_update = reshape_outputs(rates_output.gsaer_update,trim_empty);
rates_output_by_sat.gslink_rates_update = reshape_outputs(rates_output.gslink_rates_update,trim_empty);

trim_empty = true;
viz_output_by_sat.t_o = reshape_outputs(viz_output.t_o,trim_empty);
viz_output_by_sat.o_locations = reshape_outputs(viz_output.o_locations,trim_empty);
viz_output_by_sat.t_d = reshape_outputs(viz_output.t_d,trim_empty);
viz_output_by_sat.d_part = reshape_outputs(viz_output.d_part,trim_empty);
viz_output_by_sat.t_x = reshape_outputs(viz_output.t_x,trim_empty);
viz_output_by_sat.x_part = reshape_outputs(viz_output.x_part,trim_empty);

trim_empty = false;
viz_output_by_sat.dlnk_rate_history = reshape_outputs(viz_output.dlnk_rate_history,trim_empty);
viz_output_by_sat.xlnk_rate_history = reshape_outputs(viz_output.xlnk_rate_history,trim_empty);
