# Perform the job

import json
import subprocess


def call_job(params_json_file="params_template.json"):
    #output_filename = params_json_file.replace(".json", ".out")
    #subprocess.check_call("python solve_full_dynamics.py -i %s -p false -s %s" % (params_json_file, output_filename), shell=True)
    output_filename = params_json_file.replace(".json", "_fulldm.out")
    subprocess.check_call("python solve_full_dynamics_fulldm.py -i %s -p false -s %s" % (params_json_file, output_filename), shell=True)
    output_filename = params_json_file.replace(".json", "_matrixproduct.out")
    subprocess.check_call("python solve_full_dynamics_dm_matrixproduct.py -i %s -p false -s %s" % (params_json_file, output_filename), shell=True)


#call_job(params_json_file="params_N_1.json")
#call_job(params_json_file="params_N_2.json")
#call_job(params_json_file="params_N_3.json")
#call_job(params_json_file="params_N_4.json")
#call_job(params_json_file="params_N_6.json")
call_job(params_json_file="params_N_3_angle_2.5.json")
