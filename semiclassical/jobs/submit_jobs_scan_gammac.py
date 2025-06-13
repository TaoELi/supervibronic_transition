# Perform the job

import json
import subprocess


params_json_file = "params_template.json"

# Scan over the molecular number
gammac_lst = [1e-6, 3e-6, 1e-5, 2e-5, 5e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
for gammac in gammac_lst:
    print("Scanning gamma_c = ", gammac)
    # Load the parameters from the json file
    with open(params_json_file, "r") as handle:
        params = json.load(handle)
    params["gamma_c"] = gammac
    # Write the parameters to the json file
    params_json_file_tmp = "params_gammac_%.0e.json" % gammac
    print("Writing the parameters to the json file", params_json_file_tmp)
    with open(params_json_file_tmp, "w") as outfile:
        json_object = json.dumps(params, indent=4)
        outfile.write(json_object)
    # After writing the json file, submit the job
    output_filename = params_json_file_tmp.replace(".json", ".out")
    subprocess.check_call("python solve_full_dynamics.py -i %s -p false -s %s" % (params_json_file_tmp, output_filename), shell=True)
