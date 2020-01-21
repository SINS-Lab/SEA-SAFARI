import os           # path/file searching, etc
import subprocess   # Subprocesses for each run
import time         # Delay between spawning of runs
import platform     # Linux vs Windows checks

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
            file = os.path.join(directory, filename).replace('.input','')
            saf_file = open('safari.input', 'w')
            saf_file.write(file)
            saf_file.close()
            print('Running: '+filename)
            process = subprocess.Popen(command, shell=True)
            processes.append(process)
            #Wait a second for it to start running
            time.sleep(1)
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
#    run()
    txtnum = input("Number of Threads? ")
    rundir = input("Run Directory? ")
    threads = int(txtnum)
    if rundir!='.':
        rundir = os.path.join('.',rundir)
    processes = []
    run(threads, processes, directory=rundir)
    # Wait for all remaining processes to finish before closing.
    for p in processes:
        p.wait()