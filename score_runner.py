#########################
# Author: Arne Roeters 	#
# Bioinformatics group 	#
# Wageningen			#
#########################
import sys, os, argparse
import subprocess
import glob
import numpy as np
import threading
import time

def get_cmd_arguments():
	"""
	Retrieves the arguments from the command line that are input for the algorithm
	"""
	parser = argparse.ArgumentParser(prog="score_runner.py")
	parser.add_argument("-b", "--compound_file", help="The bgc to compound file, with complete path to the file",required=True)
	parser.add_argument("-g", "--group_file", help="The compound to class file, with complete path to the file",required=True)
	parser.add_argument("-n", "--network_file_folder", help="The BGC network file folder or a single file",required=True)
	parser.add_argument("-o", "--output_dir", help="The directory to put the resulting score file in.", required=False, default="./")
	parser.add_argument("-f", "--filtered_network", help="If used, the program will also output a filtered network file based on the affinity propagation.", action='store_true')
	parser.add_argument("-jw", help="The weight for the jaccard index (multiple weights seperate by , should be same length as aiw,ddsw and anchor_boost)", required=False, default="0.33")
	parser.add_argument("-aiw", help="The weight for the adjacency index, (multiple weights seperate by , should be same length as jw,ddsw and anchor_boost)", required=False, default="0.33")
	parser.add_argument("-dssw", help="The weight for the domain duplication score, (multiple weights seperate by , should be same length as jw,aiw and anchor_boost)", required=False, default="0.34")
	parser.add_argument("-anchor_boost", help="The weight for the anchor boost, (multiple weights seperate by , should be same length as jw, aiw and ddsw)", required=False, default="1.0")
	parser.add_argument("-new", help="Use the new (T) or the old (F, default) affinity propagation", required=False, default="F")
	parser.add_argument("-AP", help="if you want to do the AP add the flag", required=False, action='store_true')
	parser.add_argument("-FS", help="if you want to do the network scoring add the flag", required=False, action='store_true')
	
	args = parser.parse_args()
	
	if os.path.isdir(args.network_file_folder) and os.path.isdir(args.output_dir) and os.path.isfile(args.compound_file) and os.path.isfile(args.group_file):
		pass
	elif os.path.isfile(args.network_file_folder) and os.path.isdir(args.output_dir) and os.path.isfile(args.compound_file) and os.path.isfile(args.group_file):
		pass
	else:
		sys.exit("check all files and folders you specified for typos")
		
	try:
		new_ap = args.new.upper()
		if new_ap != "T" and new_ap != "F":
			sys.exit("please enter a T or F for the new or old AP")
	except:
		sys.exit("please enter a T or F for the new or old AP\nEntered value: " + args.new)
		
	jw = args.jw.split(",")
	aiw = args.aiw.split(",")
	dssw = args.dssw.split(",")
	anchor = args.anchor_boost.split(",")
	n = len(jw)
	if all(len(x) == n for x in [aiw,dssw,anchor]):
		return args.compound_file, args.group_file, args.network_file_folder, args.output_dir, args.filtered_network, jw, aiw, dssw, anchor, new_ap, args.AP, args.FS
	else:
		sys.exit("Check the length of the lists of all you weights")


		
def main():
	cpf, gf, nff, out, ff, jw, aiw, ddsw, anchor_boost, new_ap, do_AP, do_FS = get_cmd_arguments()
	if os.path.isdir(nff):
		network_files = glob.glob(nff + "/*.network")
		network_files.extend(glob.glob(nff + "/*.network.bz2"))
	elif os.path.isfile(nff):
		network_files = list()
		network_files.append(nff)
	if len(network_files) == 0:
		sys.exit("no files found")
	print("\nFound files:\n" + str(network_files) + "\n")
	
	threads = [] # will store all threads that have to be finished.

	"""
	This is used to do a large number of different weights all at once in a multithreaded way		
	"""
	# create a list with all processes in it
	thread_list = list()
	for Jw in np.arange(0.0, 1.05, 0.05):
		for AIw in np.arange(0.0, 1.05, 0.05):
			for DSSw in np.arange(0.0, 1.05, 0.05):
				for anchor_boost in np.arange(1.0, 5.0):
					# would for some reason sometimes add 0.0000000000000001
					Jw = round(Jw, 2)
					AIw = round(AIw, 2)
					DSSw = round(DSSw, 2)
					anchor_boost = round(anchor_boost, 2)
					if (Jw + AIw + DSSw) == 1.0:
						for nf in network_files:
							if new_ap == "F":
								propagation_file = "{}/affinity_propagation_{}_Jw{}_AIw{}_DSSw{}_AB{}.csv".format(out, nf.split("/")[-1].replace(".network", "").replace(".bz2", ""), Jw, AIw, DSSw, anchor_boost)
							else:
								propagation_file = "{}/new_affinity_propagation_{}_Jw{}_AIw{}_DSSw{}_AB{}.csv".format(out, nf.split("/")[-1].replace(".network", "").replace(".bz2", ""), Jw, AIw, DSSw, anchor_boost)
							if new_ap == "F":
								command_ap = "python Affinity_Popagation_new_weights.py -n {} -j {} -a {} -d {} -b {} -o {}/".format(nf, Jw, AIw, DSSw, anchor_boost, out)
							else:
								command_ap = "python Affinity_Popagation_new_weights.py -n {} -j {} -a {} -d {} -b {} -o {}/ -p".format(nf, Jw, AIw, DSSw, anchor_boost, out)
							if ff:
								command_nfs = "python NFS.py -n {} -a {} -o {} -g {} -b {} -f".format(nf, propagation_file, out, gf, cpf)
							else:
								command_nfs = "python NFS.py -n {} -a {} -o {} -g {} -b {} ".format(nf, propagation_file, out, gf, cpf)
							t = MyThread(command_ap, command_nfs, propagation_file, do_AP, do_FS)
							thread_list.append(t)
	
	# Run all processes that are made						
	thread_nr = 0					
	total_threads = len(thread_list)
	while True:
		if threading.activeCount() < 12 and thread_nr < total_threads:
			t = thread_list[thread_nr]
			print("starting " + t.getID())
			threads.append(t)
			t.start()
			thread_nr += 1
	# Flag = 1 means thread is still alive, 0 means no threads alive anymore
	flag = 1
	while (flag):
		time.sleep(0.5)
		flag = isThreadAlive(threads)
	
	
class MyThread(threading.Thread):
	"""
	A class to use for multithreading, extends the already existing threading.Thread
	"""
	def __init__(self, command_ap, command_nfs, ap, do_AP, do_FS):
		"""
		Initiates the class with the given arguments
		"""
		threading.Thread.__init__(self)
		self.command_ap = command_ap
		self.command_nfs = command_nfs
		self.ThreadID = ap
		self.do_AP = do_AP
		self.do_FS = do_FS
		
	def run(self):
		"""
		When Thread.start() is called on a Thread object, this will be called by the start() method in turn.
		"""
		if self.do_AP:
			subprocess.call(self.command_ap, shell=True)
		if self.do_FS:
			subprocess.call(self.command_nfs, shell=True)
		print("Finished scoring of: " + self.ThreadID)
	
	def getID(self):
		return self.ThreadID

def isThreadAlive(threads):	
	"""
	Check if any threads are still alive
	"""
	for t in threads:
		if t.isAlive():
			return 1
	return 0	

if __name__ == "__main__":
	main()
	print("########################\n# Finished the scoring #\n########################")
