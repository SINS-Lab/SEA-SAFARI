#!/usr/bin/env python3

import os           # path/file searching, etc
import subprocess   # Subprocesses for each run
import time         # Delay between spawning of runs
import platform     # Linux vs Windows checks
import argparse     # parses command line arguments

def run(threads, processes, directory='.', recursive=True):
    
    command = 'Sea-Safari.exe'
    if platform.system() == 'Linux':
        command = './Sea-Safari'
    
    for filename in os.listdir(directory):
        truefile = os.path.join(directory, filename)
        if recursive and os.path.isdir(truefile):
            run(threads, processes, truefile)
            continue
        if filename.endswith(".input") and not filename==('safari.input'):
            arg = os.path.join(directory, filename).replace('.input','')
            run_command = command + " -i " + arg
            print('Running: '+filename)
            process = subprocess.Popen(run_command, shell=True)
            processes.append(process)
            #Increment our counter of running processes
            number = len(processes) + 1
            # If hit number to run, wait for a run to finish, before continuing.
            while number > threads:
                for p in processes:
                    p.poll()
                    if p.returncode != None:
                        processes.remove(p)
                        break
                number = len(processes) + 1
                time.sleep(1)
    return
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--count", help="Number of processes")
    parser.add_argument("-d", "--dir", help="Base directory to process")
    args = parser.parse_args()
    txtnum = args.count
    rundir = args.dir
    if txtnum is None:
        txtnum = input("Number of Threads? ")
    if rundir is None:
        rundir = input("Run Directory? ")
    threads = int(txtnum)
    if rundir!='.':
        rundir = os.path.join('.',rundir)
    processes = []
    run(threads, processes, directory=rundir)
    # Wait for all remaining processes to finish before closing.
    for p in processes:
        p.wait()