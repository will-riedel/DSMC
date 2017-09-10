# test
# Will Riedel
# 2017

import numpy as np
from datetime import datetime
import sys
from subprocess import call
import time

def main(argv):

	



	t_alloted = int(5)
	mem_alloted = int(1000)


	# command_str = ("python DSMC_input.py "+
	command_str = ("./Main")
	print command_str

	sbatch_str = ("-p normal -t "+
					str(t_alloted)+
					" --mem="+str(mem_alloted)+
					" --nodes=1 --job-name=DSMC --output=out/DSMC.%j.out --error=err/DSMC.%j.err --mail-type=FAIL --mail-user=wriedel@stanford.edu"
					)
	print sbatch_str

	output_str = 'sbatch '+sbatch_str+' --wrap="'+command_str+'"'
	call("cd /home/wriedel/DSMC_FORTRAN",shell=True)
	call(output_str,shell=True)
	time.sleep(1)





if __name__ == "__main__":
	main(sys.argv)
















print "----- done -----"

