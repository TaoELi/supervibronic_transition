# Perform the job

import json
import subprocess


def call_job(params_json_file="params_template.json"):
    output_filename = params_json_file.replace(".json", ".out")
    subprocess.check_call("python solve_full_dynamics.py -i %s -p false -s %s" % (params_json_file, output_filename), shell=True)


#call_job(params_json_file="params_template.json")
#call_job(params_json_file="params_decoupled_e.json")
#call_job(params_json_file="params_decoupled_e_withfield.json")
call_job(params_json_file="params_decoupled_eph.json")
#call_job(params_json_file="params_decoupled_phv.json")
#call_job(params_json_file="params_decoupled_phv_noloss.json")
#call_job(params_json_file="params_decoupled_phv_nocavloss.json")
#call_job(params_json_file="params_decoupled_phv_nodarkloss.json")
