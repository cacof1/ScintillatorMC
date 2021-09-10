import os
from random import seed
from random import randint
for i in range(0,100):
    seed = randint(0,1000000)
    f = open('command.sh','w')
    f.write('#!/bin/bash -l \n')
    f.write('#$ -l h_rt=12:00:0 \n')
    f.write('#$ -l mem=5G \n')
    f.write('#S -wd /home/cacfek/ScintillatorMC \n')
    f.write('source /home/cacfek/software/root/bin/thisroot.sh \n')
    f.write('cd /home/cacfek/ScintillatorMC \n')
    f.write('/home/cacfek/ScintillatorMC/bin/ScintillatorMC pCT_config.txt ' + str(seed))
    f.close()
    os.system('qsub command.sh')
    os.system('rm command.sh')
