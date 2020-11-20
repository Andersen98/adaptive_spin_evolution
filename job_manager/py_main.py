from adaptive_io import mode_generator as mg
import subprocess
from collections import defaultdict
import linecache
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
def generate_and_submit(doRun):
    
    BATCH_N = 16;

    #Run Info
    prefix_code = "/home/ethan/curr_proj/run_output/"
    executable = "/home/ethan/curr_proj/adaptive_spin"
    run_id=1428
    
    #Run Info (system count)
    num_modes = 1000
    max_occupation = 15

    #Energy Paths
    prefix_energy = "/home/ethan/curr_proj/energy_params/"
    a_name = "atom_energies"
    m_name = "mode_energies"
    
    #Energy Params
    remake_energies = False
    w0 = 1
    g0 = 10
    v0 = 1
    cutoff = .00001
    spectral_energy = 1

    #Time Params
    dt = .001
    tf =2000*dt
    
    #=================================================================#
    #                      Make Energies                               #

    
    #=================================================================#
    #                     Generate Input Energies                     #

    mode_gen = mg.mode_generator(w0,g0,v0,spectral_energy,num_modes,dt)    
    w,g = mode_gen.w_g_rnd_dist()

    #                      End Subsection                             #
    #=================================================================#


    
    #=================================================================#
    #                     Write Input Energies                        #

    #arguments for adaptive spin
    with open(prefix_energy+ m_name,"w") as f:
        for i in range(num_modes):
            #don't add a new line if it is the last one
            append = "\n" if i < (num_modes-1) else ""
            f.write( str(w[i]) +"\t")
            f.write( str(g[i]) +append)

    with open(prefix_energy+ a_name,"w") as f:
        f.write(str(0) + "\n")
        f.write(str(w0))

    #                      End Subsection                             #
    #=================================================================#

        
    #                         End Section                             #
    #=================================================================#


    #=================================================================#
    #                Prepare Arguments for Adaptive Spin              #
    arguments = {

        "--run_id=":run_id,
        "--output_dir=":prefix_code,
        "-a":prefix_energy+a_name,
        "-m":prefix_energy+m_name,
        "--cutoff=":cutoff,
        "--tf=":tf,
        "--dt=":dt
        }


    argList = [executable]
    for key,value in arguments.items():
        argList += [str(key)+ str(value)]
    
         
    #                         End Section                             #
    #=================================================================#

    

    #=================================================================#
    #                    Run Instance of C++ Code                     #
    
    if(doRun):
        subprocess.run(argList,check=True)
    

    #                         End Section                             #
    #=================================================================#
    outStr = ""
    for x in argList:
        outStr += x + " "

    print( outStr)

    return(arguments,w,g,w0)
 
def p_excited_strong(t,K,G,D,g0):
    
    g = np.sqrt(g0**2+(D/2)**2-(G/2)**2)
    factor = np.exp(-K*t)/2
    term1 = (G**2 +D**2)/(4*g**2)
    term2 = (1-(G**2 +D**2)/(4*g**2))*np.cos(2*g*t)
    term3 = G/g * np.sin(2 *g* t)
    result =factor*(1+ term1 +term2 +term3)
    return result
def p_excited(t,w,g,w0):
    idx = np.argmin(np.abs(g-w0))
    g_near = g[idx]
    volume = 1/(np.abs(w[1]-w[0]))
    Gamma = w0**3
    print("Gamma",Gamma)
    return np.exp(-Gamma *t)

if __name__=='__main__':

    arguments,w,g,w0 = generate_and_submit(True)

    run_file = arguments["--output_dir="]+ str( arguments["--run_id="])+".out"
    print(run_file)
    locations = []
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
        print(header)

        while(count < locations["DSTART"]-1):
            f.readline()
            count = count +1
        
        
        data_list = np.empty((1,2),float)
        for data in f:
            spl= data.split()
            run_id = spl[0]
            print(spl)
            val = np.empty((1,2),float)
            val[0]   = np.array([float(spl[1]),float(spl[2])],float)
            print(val)
            data_list = np.append(data_list,val,axis=0)

        data_list = data_list[1:,:]

        title = "Plot of Spin Populations"
        xLbl = "time ( O(1/w0))"
        yLbl = "Population"
        fig, ax = plt.subplots()
        t = data_list[:,0]
        up = data_list[:,1]
        tx = np.linspace(0,t[-1],100)
        print(data_list)
        print(p_excited(tx,w,g,w0))
        plt.semilogy(tx,p_excited(tx,w,g,w0),label="Theory:~exp(-omega_0^3 t")
        plt.semilogy(t,up,label="Code")
        plt.grid()
        plt.legend()
        ax.set(xlabel=xLbl,ylabel=yLbl,title=title)
        fig.savefig("./figs/spin_pop"+str( arguments["--run_id="])+".png")
        plt.close()
        
        
        
