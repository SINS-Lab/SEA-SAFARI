#!/usr/bin/env python3

import subprocess   # used to run safari and vmd
import argparse     # parses command line arguments
import platform     # Linux vs Windows check
import os           # os.remove is used for .vmd file
import time         # sleeps before removing file
import hashlib      # Used to check if the files matched the old hashes


# These are the expected file hashes after running a test.
# If you make a change which affects output, you need
# to update these file hashes as well, otherwise
# the automated build will fail!
MonteCarloHash = '085cfa7165309e489ee805873c1cac4a'
AdaptiveGridHash = '29868be762053fc6a1d481ba4066ec80'
ChainScatHash = 'e9c5057efd45618e460a0246ba78c5c9'



def test_montecarlo(template):
    command = template.format('test_montecarlo', 'tests/test_montecarlo')
    print("Running Test of Montecarlo\n")
    subprocess.run(command, shell=True)
    md5_found = hashlib.md5(open('./tests/test_montecarlo.spec','rb').read()).hexdigest()
    if (MonteCarloHash != md5_found):
        print("Hash Mis-match for Monte Carlo! New hash is:")
        print(md5_found)
        return 1
    else:
        print('Monte Carlo hash matched correctly!')
        return 0
    
def test_adaptive_grid(template):
    command = template.format('test_adaptive_grid', 'tests/test_adaptive_grid')
    print("\n\nRunning Test of Adaptive Grid\n")
    subprocess.run(command, shell=True)
    md5_found = hashlib.md5(open('./tests/test_adaptive_grid.spec','rb').read()).hexdigest()
    if (AdaptiveGridHash != md5_found):
        print("Hash Mis-match for Adaptive Grid! New hash is:")
        print(md5_found)
        return 2
    else:
        print('Adaptive Grid hash matched correctly!')
        return 0

def test_chainscat(template):
    command = template.format('test_chain', 'tests/test_chain')
    print("\n\nRunning Test of Chainscat\n")
    subprocess.run(command, shell=True)
    md5_found = hashlib.md5(open('./tests/test_chain.spec','rb').read()).hexdigest()
    if (ChainScatHash != md5_found):
        print("Hash Mis-match for Chain Scat! New hash is:")
        print(md5_found)
        return 4
    else:
        print('Chain Scat hash matched correctly!')
        return 0

if __name__ == '__main__':
    # Run a single shot safari for this run,
    # assuming the input file 
    # was already configured properly.
    template = 'Sea-Safari.exe -i {} -o {}'

    #Change command accordingly for linux
    if platform.system() == 'Linux':
        template = template.replace('Sea-Safari.exe', './Sea-Safari')

    if not os.path.exists('./tests'):
        os.makedirs('./tests')

    result = test_montecarlo(template)
    result = result + test_adaptive_grid(template)
    result = result + test_chainscat(template)

    if result != 0:
        print("Tests Failed, or a Hash was intentionally changed, if this was intentional, please update the expected hashes!")
        raise SystemExit('Error: Hash mis-match was found!')
    
    print("Finished running all Tests!")