import subprocess
import os
import time
import glob

def main():
	"""
	A quick build antismash runner
	"""
	files = glob.glob("/mnt/scratch/roete009/antismash_input/*")
	processes = set()
	max_processes = 12
	for name in files:
		dir_name = name.split("/")[5].split(".")[0]
		subprocess.call("mkdir /mnt/scratch/roete009/antismash_output/{}".format(dir_name), shell=True)
		command = "antismash  --pfamdir /mnt/scratch/roete009/pfam/ --smcogs --borderpredict --outputfolder /mnt/scratch/roete009/antismash_output/{} --logfile /mnt/scratch/roete009/antismash_output/{}/log_file.txt -c 12 {}".format(dir_name, dir_name, name)
		command = command.split()
		processes.add(subprocess.Popen(command))
		if len(processes) >= max_processes:
			os.wait()
			processes.difference_update([p for p in processes if p.poll() is not None])
	
if __name__ == "__main__":
	main()




