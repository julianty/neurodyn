#!/share/apps/opt/python/2.7.9/bin/python2
#$ -S /share/apps/opt/python/2.7.9/bin/python2
#$ -V
#$ -cwd
#$ -j y
#$ -M prozdeba@physics.ucsd.edu
#$ -o ./output
#$ -e ./error
#$ -q batch.q

import os


Ninit = 10

SGE_TASK_ID = int(os.getenv("SGE_TASK_ID", 0))

initID = int(SGE_TASK_ID - 1) % Ninit + 1
adolcID = SGE_TASK_ID % 2000

print("initID = %d"%(initID,))
print("SGE_TASK_ID = %d"%(SGE_TASK_ID,))

print(os.system("uname -n"))

#if i_M == 4:
#    print("Skipping M=1000 case.")
#else:
os.system("python2 ./varanneal/VarAnneal_Neurodyn.py  %d %d"%(initID, M[i_M], D_hidden[i_DH], adolcID))
