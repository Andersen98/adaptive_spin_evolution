import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import json

def import_run_data(json_path,suffix):
    locations = dict()
    conf = dict()
    
    with open(json_path,"r") as f:
        json_string = f.read(-1)
        conf = json.loads(json_string)

    run_paths = conf["run_info"]["system_paths"]
    run_id = conf["run_info"]["run_id"]
    run_file = run_paths["code_output_dir"] + str(run_id) + "."+suffix
    with open(run_file,"r") as f:
        count = 0
        data_header = f.readline()
        count = count +1
        data_header = data_header.split()
        locations = {el.split("~")[0]:int(el.split("~")[1]) for el in data_header }

        header = ""
        while(count < locations["PSTART"]-1):
            f.readline()
            count = count +1
        header = f.readline()
        count = count +1

        while(count < locations["DSTART"]-1):
            f.readline()
            count = count +1
        
        
        data_list = np.empty((1,3),float)
        for data in f:
            spl= data.split()
            run_id = spl[0]
            val = np.empty((1,3),float)
            val[0]   = np.array([float(spl[1]),float(spl[2]),float(spl[3])],float)
            data_list = np.append(data_list,val,axis=0)

        data_list = data_list[1:,:]
        return data_list,conf

def plot_run(json_path,dunp_data=False):
    data_list,conf = import_run_data(json_path,"full_spin_pop")
    title = "Plot of Spin Populations"
    xLbl = "time (hbar/eV ~6.5 *10^-16 sec)"
    yLbl = "Population"
    fig, ax = plt.subplots()
    t = data_list[:,0]
    up = data_list[:,1]
    #print(data_list)
    
    
    plt.plot(t,up,label="Code")
    plt.grid()
    plt.legend()
    ax.set(xlabel=xLbl,ylabel=yLbl,title=title)
    run_paths = conf["run_info"]["system_paths"]
    fig.savefig(run_paths["code_output_dir"]+"figs/"+str(conf["run_info"]["run_id"])+"full_spin_pop.png")
    plt.close()

    data_list,conf = import_run_data(json_path,"cavity_emitter_pair")
    title = "Plot of Cavity Mode Probability (1 quanta)"
    xLbl = "time (hbar/eV ~6.5 *10^-16 sec)"
    yLbl = "Probability"
    fig, ax = plt.subplots()
    t = data_list[:,0]
    #print(data_list)
    emitter = data_list[:,1]
    cavity = data_list[:,2]
    
    plt.plot(t,emitter,label="Emitter Excited, Cavity Ground")
    plt.plot(t,cavity,label="Cavity Excited,Emitter Ground")
    plt.grid()
    plt.legend()
    ax.set(xlabel=xLbl,ylabel=yLbl,title=title)
    run_paths = conf["run_info"]["system_paths"]
    fig.savefig(run_paths["code_output_dir"]+"figs/"+str(conf["run_info"]["run_id"])+"cavity_emitter_pair.png")
    plt.close()
