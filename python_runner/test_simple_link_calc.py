import json
from simple_link_calc import pack_links_data

test_case = 'zhou'
abs_path = r'C:\Users\User\circinusGit\SPRINT\source\circinus_orbit_link_public\python_runner\example_output'
data_rates_filename = abs_path+'\\data_rates_output_%s.json'%test_case
orbit_prop_filename = abs_path +'\\orbit_prop_data_%s.json'%test_case

with open(data_rates_filename,'r') as f:
        data_rates_example = json.load(f)

with open(orbit_prop_filename,'r') as f:
        orbit_prop_data = json.load(f)
# Only testing data_rates output
""" sat_link_history_filename = rel_path+'sat_link_history_%s.json'%test_case
with open(sat_link_history_filename,'r') as f:
        sat_link_history = json.load(f) """


dlnk_rate_Mbps = 20
xlnk_rate_Mbps = 10

# MAX duration for crosslinks and downlinks (TODO: enforce for downlinks, only enforce Xlnks now)
xlnk_max_len_s = 200

example_data_rates = data_rates_example['accesses_data_rates']

test_data_rates = pack_links_data(orbit_prop_data,dlnk_rate_Mbps,xlnk_rate_Mbps,xlnk_max_len_s)

print('all data loaded, ready for testing')

# TODO: write some tests, I just did some basic inspections for now 