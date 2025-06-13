# Perform the job

import json
import subprocess


params_json_file = "params_template.json"

# Scan over the incoming electric field amplitude
E0_ext_lst = [1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2, 2e-2, 4e-2]
for E0_ext in E0_ext_lst:
    print("Scanning E0_ext = ", E0_ext)
    # Load the parameters from the json file
    with open(params_json_file, "r") as handle:
        params = json.load(handle)
    params["E0_ext"] = E0_ext
    # Write the parameters to the json file
    params_json_file_tmp = "params_E0_%.0e.json" % E0_ext
    print("Writing the parameters to the json file", params_json_file_tmp)
    with open(params_json_file_tmp, "w") as outfile:
        json_object = json.dumps(params, indent=4)
        outfile.write(json_object)
    # After writing the json file, submit the job
    output_filename = params_json_file_tmp.replace(".json", ".out")
    subprocess.check_call("python solve_full_dynamics.py -i %s -p false -s %s" % (params_json_file_tmp, output_filename), shell=True)
