# Perform the job

import json
import subprocess

import numpy as np

params_json_file = "params_template.json"

angle_lst = [np.pi]
for angle in angle_lst:
    print("Scanning angle threshold = ", angle * 180 / np.pi, " degrees")
    # Load the parameters from the json file
    with open(params_json_file, "r") as handle:
        params = json.load(handle)
    params["v_dipole_angle_thres"] = angle
    # Write the parameters to the json file
    params_json_file_tmp = "params_angle_%.3f.json" % angle
    print("Writing the parameters to the json file", params_json_file_tmp)
    with open(params_json_file_tmp, "w") as outfile:
        json_object = json.dumps(params, indent=4)
        outfile.write(json_object)
    # After writing the json file, submit the job
    output_filename = params_json_file_tmp.replace(".json", ".out")
    subprocess.check_call("python solve_full_dynamics.py -i %s -p false -s %s" % (params_json_file_tmp, output_filename), shell=True)
    output_filename = params_json_file_tmp.replace(".json", "_vangle.out")
    subprocess.check_call("python solve_full_dynamics_vib_angle.py -i %s -p false -s %s" % (params_json_file_tmp, output_filename), shell=True)
