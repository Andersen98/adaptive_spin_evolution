from adaptive_io import mode_generator as mg
import subprocess


if __name__=='__main__':

    
    
    BATCH_N = 16;

    #Run Info
    prefix_code = "/home/ethan/curr_proj/run_output/"
    executable = "/home/ethan/curr_proj/adaptive_spin"
    run_id=123
    #Run Info (system count)
    num_modes = 2
    max_occupation = 15

    #Energy Paths
    prefix_energy = "/home/ethan/curr_proj/energy_params/"
    a_name = "atom_energies"
    m_name = "mode_energies"
    
    #Energy Params
    remake_energies = False
    w0 = 1.5
    g0= 0
    v0 = 1.5
    cutoff =0;
    spectral_energy = .01

    #Time Params
    dt = .001
    tf =20*dt
    
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
        "--dt=":dt,
        "-v":""
        }


    argList = [executable]
    for key,value in arguments.items():
        argList += [str(key)+ str(value)]
    
         
    #                         End Section                             #
    #=================================================================#

    

    #=================================================================#
    #                    Run Instance of C++ Code                     #
    
    
    subprocess.run(argList,check=True)
    

    #                         End Section                             #
    #=================================================================#
    outStr = ""
    for x in argList:
        outStr += x + " "

    print( outStr)
    
