# Perform the job

import json
import subprocess


params_json_file = "params_template.json"

# Scan over the cavity mode frequency
omega_e_lst = [0.015, 0.020, 0.025, 0.035, 0.040, 0.050, 0.060, 0.070, 0.080, 0.100, 0.12, 0.14, 0.16, 0.20]

for omega_e in omega_e_lst:
    print("Scanning omega_e = ", omega_e)
    # Load the parameters from the json file
    with open(params_json_file, "r") as handle:
        params = json.load(handle)
    params["omega_e"] = omega_e
    params["mue"] = 0.0
    print("the permanent dipole moment is set to zero for the electronic excited state")
    # Write the parameters to the json file
    params_json_file_tmp = "params_omegae_%.3f_mue0.json" % omega_e
    print("Writing the parameters to the json file", params_json_file_tmp)
    with open(params_json_file_tmp, "w") as outfile:
        json_object = json.dumps(params, indent=4)
        outfile.write(json_object)
    # After writing the json file, submit the job
    output_filename = params_json_file_tmp.replace(".json", ".out")
    subprocess.check_call("python solve_full_dynamics.py -i %s -p false -s %s" % (params_json_file_tmp, output_filename), shell=True)
