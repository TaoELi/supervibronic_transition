# Perform the job

import json
import subprocess


params_json_file = "params_template.json"

# Scan over the molecular number
#N_lst = [1e8, 3e8, 1e9, 3e9, 1e10, 3e10, 2e8, 4e8, 5e8, 6e8, 8e8]
N_lst = [2e10, 3e10, 4e10, 6e10]
for N in N_lst:
    print("Scanning N = ", N)
    # Load the parameters from the json file
    with open(params_json_file, "r") as handle:
        params = json.load(handle)
    params["Nv"] = N
    # Write the parameters to the json file
    params_json_file_tmp = "params_Nv_%.0e.json" % N
    print("Writing the parameters to the json file", params_json_file_tmp)
    with open(params_json_file_tmp, "w") as outfile:
        json_object = json.dumps(params, indent=4)
        outfile.write(json_object)
    # After writing the json file, submit the job
    output_filename = params_json_file_tmp.replace(".json", ".out")
    subprocess.check_call("python solve_full_dynamics.py -i %s -p false -s %s" % (params_json_file_tmp, output_filename), shell=True)
