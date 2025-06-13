# Perform the job

import json
import subprocess


params_json_file = "params_template.json"

# Scan over the cavity mode frequency
omega_c_lst = [0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.011, 0.012, 0.013, 0.014, 0.015]

for omega_c in omega_c_lst:
    print("Scanning omega_c = ", omega_c)
    # Load the parameters from the json file
    with open(params_json_file, "r") as handle:
        params = json.load(handle)
    params["omega_c"] = omega_c
    # Write the parameters to the json file
    params_json_file_tmp = "params_omegac_%.3f.json" % omega_c
    print("Writing the parameters to the json file", params_json_file_tmp)
    with open(params_json_file_tmp, "w") as outfile:
        json_object = json.dumps(params, indent=4)
        outfile.write(json_object)
    # After writing the json file, submit the job
    output_filename = params_json_file_tmp.replace(".json", ".out")
    subprocess.check_call("python solve_full_dynamics.py -i %s -p false -s %s" % (params_json_file_tmp, output_filename), shell=True)
