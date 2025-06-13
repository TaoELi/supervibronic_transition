# Perform the job

import json
import subprocess


params_json_file = "params_template.json"

# Scan over the excited-state permanent dipole moment
mue_lst = [0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]

for mue in mue_lst:
    print("Scanning mue = ", mue)
    # Load the parameters from the json file
    with open(params_json_file, "r") as handle:
        params = json.load(handle)
    params["mue"] = mue
    # Write the parameters to the json file
    params_json_file_tmp = "params_mue_%.1f.json" % mue
    print("Writing the parameters to the json file", params_json_file_tmp)
    with open(params_json_file_tmp, "w") as outfile:
        json_object = json.dumps(params, indent=4)
        outfile.write(json_object)
    # After writing the json file, submit the job
    output_filename = params_json_file_tmp.replace(".json", ".out")
    subprocess.check_call("python solve_full_dynamics.py -i %s -p false -s %s" % (params_json_file_tmp, output_filename), shell=True)
