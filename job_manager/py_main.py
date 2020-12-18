from adaptive_io import generate_config as gc
from collections import defaultdict
from adaptive_io import make_plot as mp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
import os
if __name__=='__main__':

    

    params = {
        "run_info":{"run_id":1429,"num_modes":2,"system_paths":
                    {"code_output_dir":"/home/ethan/code/run_output/1429/",
                     "excecutable":"/home/ethan/code/adaptive_spin",
                     }
                    },
        "time_params":{"dt":.01,"tf":4},

    "initial_state":
        [{"re":.5,"im":0,"spin":True,"idx":0,"n":0},
         {"re":.5,"im":0,"spin":False,"idx":0,"n":1}],

        "energy_info":{
            "params":{
            "cutoff":1,
                "w0":3,
                "g0":.2,
                "v0":3,
                "energy_spectral_density":.001
            },
          
        }
    }
    if not os.path.exists(params["run_info"]["system_paths"]["code_output_dir"]):
        os.makedirs(params["run_info"]["system_paths"]["code_output_dir"])
    if not os.path.exists(params["run_info"]["system_paths"]["code_output_dir"]+"figs/"):
        os.makedirs(params["run_info"]["system_paths"]["code_output_dir"]+"figs/")
        
    arg_list,conf,json_path = gc.generate_json(params)

   
    
    
    print(subprocess.run(arg_list,check=True))

    print(json_path)
    mp.plot_run(json_path)     
        
        
