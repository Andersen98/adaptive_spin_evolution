import subprocess
import json

import random
import numpy as np
from numpy.random import default_rng


class mode_generator:

    def __init__(self,w0_,g0_,v0_,spectral_energy_,num_modes_,dt_):

        self.w0=w0_
        self.g0=g0_
        self.v0=v0_
        self.spectral_energy=spectral_energy_
        self.num_modes=num_modes_
        self.rng=default_rng()
        self.num_modes=num_modes_
        self.fast_time = 2*np.pi/(10*w0_)
        self.long_time = 5*2*np.pi *np.sqrt(1.0/num_modes_)/spectral_energy_
    def w_g_rnd_dist(self):
        N = self.num_modes
        gamma = self.spectral_energy
        v0 = self.v0
        g0 = self.g0
        g = np.array([])
        w = np.array([])
        if(self.num_modes > 1):
            w = np.linspace(self.fast_time,self.long_time,num=(N-1))
            w = 1/w;
            g = self.rng.random(N-1)
            norm = np.sum(g**2,axis=0)
            g = gamma*g*np.sqrt(1/norm)
        return w,g




    

def generate_json(params,print_command=True,dump_template=False):
    result = None,None,None #conf,json_path,json_string
    if(dump_template):
        example_dict = {
            "run_info":{"run_id":1429,"num_modes":2,"system_paths":
                        {"code_output_dir":"/home/ethan/curr_proj/run_output/",
                         "excecutable":"/home/ethan/curr_proj/adaptive_spin"}
                        },
            
            "time_params":{"dt":.001,"tf":.003},

            "initial_state":
                [{"re":.9,"im":0,"spin":True,"idx":0,"n":1},
                 {"re":.1,"im":0,"spin":False,"idx":0,"n":1}],

            "energy_info":{
                "params":{
                    "cutoff":.00001,"w0":2,"g0":.2,"v0":2,"energy_spectral_density":1},
                "energies":{
                    "emitter":{"up":1,"down":-1},
                    "modes":[ {"w":1,"g":.2},{"w":3, "g":.0003}]
                }
            }
        }
        print(json.dumps(example_dict,indent=4))
    else:

        #Run Info
        num_modes = params["run_info"]["num_modes"]

        #energy info
        w0 = params["energy_info"]["params"]["w0"]
        v0 = params["energy_info"]["params"]["v0"]
        g0 = params["energy_info"]["params"]["g0"]
        energy_spectral_density = params["energy_info"]["params"]["energy_spectral_density"]


        #time params
        dt = params["time_params"]["dt"]
        #=================================================================#
        #                      Make Energies                               #


        energies = {"emitter":{"up":w0/2,"down":-w0/2},
                    "modes":[dict() for x in range(num_modes)]}
        
        mode_gen = mode_generator(w0,g0,v0,energy_spectral_density,num_modes,dt)    
        w,g = mode_gen.w_g_rnd_dist()
        w = np.insert(w,0,w0)
        g = np.insert(g,0,g0)
        energies["modes"] = [ {"w":wg[0],"g":wg[1]} for wg in zip(w,g)]
        params["energy_info"]["energies"] = energies

        
        #                      End Subsection                             #
        #=================================================================#
        
        
        
        #=================================================================#
        #                     Write Input Energies                        #
        
        #arguments for adaptive spin
        json_str = ""
        run_id = params["run_info"]["run_id"]
        json_out = params["run_info"]["system_paths"]["code_output_dir"] + str(run_id)+".json"
        
        with open(json_out,"w") as f:
            json_str = json.dumps(params,indent=4)
            f.write(json_str)
            print(json_str)
        print(json_out)
    
            
            
        #                      End Subsection                             #
        #=================================================================#

        
        #                         End Section                             #
        #=================================================================#

            
        #=================================================================#
        #                Prepare Arguments for Adaptive Spin              #
        arguments = {
                
            "--json_file=":json_out,
           
        }


        executable = params["run_info"]["system_paths"]["excecutable"]
        argList = [executable]

        for key,value in arguments.items():
            argList += [str(key)+ str(value)]
    
        #                         End Section                             #
        #=================================================================#

        outStr = ""
        for x in argList:
            outStr += x + " "
        if(print_command):
            print(outStr)


        #=================================================================#
        #                    Run Instance of C++ Code                     #
    
        #subprocess.run(argList,check=True)
    

        #                         End Section                             #
        #=================================================================#
    
        result = (params,json_out,json_str)

    return result
   
if __name__=="__main__":
    argList,params = generate("json_configs/test.json",dump_template=True)

