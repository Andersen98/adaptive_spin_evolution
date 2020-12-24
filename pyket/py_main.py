from collections import defaultdict
import pytools 
import pyket as pk
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import os
if __name__=='__main__':

    

    params = {
        "run_info":{"run_id":1429,"num_modes":pk.num_modes(),"system_paths":
                    {"code_output_dir":"/home/ethan/code/run_output/1429/",
                     "excecutable":"/home/ethan/code/adaptive_spin",
                     }
                    },
        "time_params":{"dt":.01,"tf":70},

        "initial_state":
        [{"re":1,"im":0,"spin":True,"idx":0,"n":0},
         {"re":0,"im":0,"spin":False,"idx":0,"n":1}],
        
        "energy_info":{
            "params":{
            "cutoff":.0004054,
                "w0":3,
                "g0":2,
                "v0":3,
                "energy_spectral_density":.01,
                
            },
          
        }
    }


        

    
    if not os.path.exists(params["run_info"]["system_paths"]["code_output_dir"]):
        os.makedirs(params["run_info"]["system_paths"]["code_output_dir"])
    if not os.path.exists(params["run_info"]["system_paths"]["code_output_dir"]+"figs/"):
        os.makedirs(params["run_info"]["system_paths"]["code_output_dir"]+"figs/")
        
    pytools.make_hamiltonian(params)

    def my_spin_pop(state_vec):
        up = 0
        down = 0
        for k in state_vec:
            if(k.get_spin()):
                up += abs(k.get_amp())**2
            else:
                down += abs(k.get_amp())**2

        return up,down
    def zero_hist(hist,n):
       
        for i in range( n-1):
            if(hist[-(1+i)] != hist[-(2+i)]):
                return False
        return True
    
    dt = .01
    h = pk.H(json_str)
    n_thresh = 4
    n_hist =[]
    n_hist.append(len(h.get_state_vector()))
    
    epsilon = .0004054
    for i in range(1000):
        
        up, down = h.get_spin_pop()
        print("C++     "+str(dt*i) + "  " + str(up) + "  " +str(down))
        h.run_step(dt)
        n_hist.append(len(h.get_state_vector()))

        if len(n_hist)> n_thresh and  zero_hist(n_hist,n_thresh):
            epsilon = epsilon*.01
            h.set_epsilon(epsilon)
            print("setting epsilon to " + str(epsilon))
            print("current state size is " + str(n_hist[-1]))
    
    #print(json_path)
    #mp.plot_run(json_path)     
        
        
