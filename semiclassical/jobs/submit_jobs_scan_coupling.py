# Perform the job

import json
import subprocess


params_json_file = "params_template.json"

# Scan over the coupling strength
#lambda_c_lst = [4e-7, 6e-7, 8e-7, 1e-6, 2e-6, 3e-6, 4e-6, 1e-8, 3e-8, 6e-8, 1e-7, 2e-7]
lambda_c_lst = [0.0]
for lambda_c in lambda_c_lst:
    print("Scanning lambda_c = ", lambda_c)
    # Load the parameters from the json file
    with open(params_json_file, "r") as handle:
        params = json.load(handle)
    params["lambda_c"] = lambda_c
    # Write the parameters to the json file
    params_json_file_tmp = "params_lambdac_%.0e.json" % lambda_c
    print("Writing the parameters to the json file", params_json_file_tmp)
    with open(params_json_file_tmp, "w") as outfile:
        json_object = json.dumps(params, indent=4)
        outfile.write(json_object)
    # After writing the json file, submit the job
    output_filename = params_json_file_tmp.replace(".json", ".out")
    subprocess.check_call("python solve_full_dynamics.py -i %s -p false -s %s" % (params_json_file_tmp, output_filename), shell=True)
