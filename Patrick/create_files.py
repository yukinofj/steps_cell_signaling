import os

runstring = "import os \n"

n = 4
R = 1.5
T0 = -100.0
P = 30.0
L = 20.0
totaltime = 100000
og_run = 2

for runnumber in range(10, 13):
    file = "n{}R{}T{}P{}L{}run{}".format(n, R, T0, P, L, og_run)
    jobscript = f"""#!/bin/bash --login
cd /home/yukinofj/BA/program/sim
/home/yukinofj/.juliaup/bin/julia "ABP_Simulation_rnd.jl" {n} {R} {T0} {P} {L} {totaltime} {runnumber} {file}
"""
    script_filename = f'jobscriptn{n}R{R}T{T0}P{P}L{L}t{totaltime}m{runnumber}.sh'
    
    with open(script_filename, "w") as text_file:
        text_file.write(jobscript)
    
    runstring += f"os.system('chmod 777 {script_filename}')\n"
    runstring += f"os.system('qsub -mem 6 -env \"PATH=$PATH;HOME=$HOME\" {script_filename}')\n\n"

with open('run_files.py', "w") as text_file:
    text_file.write(runstring)