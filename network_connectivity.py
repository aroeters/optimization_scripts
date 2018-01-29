import sys
import os
import glob
import argparse
import numpy as np
import threading
import networkx
import time

class ConnectivityCalculator:
	def __init__(self):
		"""
		Initiates the class
		"""
		self.output_dir = ""
		self.network_files = list() # network files
		self.raw_distance_components = dict()
		self.thread_list = list() # all threads
		self.threads = 30

	def getCmdArguments(self):
		"""
		Retrieves the command line arguments
		"""
		print("Retrieving command line arguments...")
		parser = argparse.ArgumentParser(prog='network_connectivity.py')
		parser.add_argument("-n", "--network_file_dir", help="The directory with the to be used network files in it", required=True)
		parser.add_argument("-o", "--output_dir", help="The directory to output the final results in", required=False, default="./")
		parser.add_argument("-c", "--threads", help="Number of threads to use (default = 30)", required=False, default=30)
		args = parser.parse_args()

		if os.path.isdir(args.output_dir):
			self.output_dir = args.output_dir
			self.open_result_file()
		else:
			sys.exit("The given output directory is not valid")

		if os.path.isdir(args.network_file_dir):
			network_files = glob.glob(args.network_file_dir + "/*.network")
			if len(network_files) != 0:
				self.network_files = network_files
			else:
				sys.exit("No network files found in the given directory")
		else:
			sys.exit("The given network directory is not valid")
		try:
			self.threads = int(args.threads)
		except:
			sys.exit("Not a valid number of threads...")
		print("Command line arguments retrieved...")

	def parse_network_files(self):
		"""
		Parses the network files and puts the different components into the dictionary
		"""
		print("Reading in all files...")
		for nw_file in self.network_files:
			self.raw_distance_components.update({nw_file:dict()})
			file_in = open(nw_file)
			file_in.readline()
			for line in file_in:
				bgc1, bgc2, raw_distance, sim_squared, Jaccard, DSS, AI, raw_dss_non_anchor, raw_dss_anchor, non_anchor_domains, anchor_domains = line.split()[0:11]
				DSS = self.calculate_DSS(1, raw_dss_non_anchor, raw_dss_anchor, non_anchor_domains, anchor_domains)
				self.raw_distance_components[nw_file].update({(bgc1, bgc2):{"J":Jaccard, "AI":AI, "DSS":DSS}})
		print("Done reading files...")

	def prepare_distance_with_weights_calculations(self):
		"""
		Calculates the distances with different weights for the Jaccard, AI and DSS
		"""
		total = 0
		for nw_file in self.raw_distance_components.keys():
			self.thread_list = list()
			print("Starting: " + nw_file)
			self.current_file = nw_file
			for Jw in np.arange(0.0, 1.05, 0.05):
				for AIw in np.arange(0.0, 1.05, 0.05):
					for DSSw in np.arange(0.0, 1.05, 0.05):
						Jw = round(Jw, 2)
						AIw = round(AIw, 2)
						DSSw = round(DSSw, 2)
						if (Jw + AIw + DSSw) == 1.0:
							self.thread_list.append(NetworkConnectingThread(Jw, AIw, DSSw, self.raw_distance_components[nw_file], nw_file, str(Jw)+str(AIw)+str(DSSw)+nw_file))
							total += 1
			print("total jobs for {} is {}".format(nw_file, str(len(self.thread_list))))
			print("Finished preparing threads...")
			self.run_calculations()
			print("Done with " + nw_file)
		print("Total number of jobs was: {}".format(total))
		

	def run_calculations(self):
		"""
		runs all the threads that are needed
		"""
		print("Starting to run threads...")
		my_current_threads = []
		total_jobs = len(self.thread_list)
		active_threads = True
		thread_nr = 0
		while total_jobs >= thread_nr or active_threads:
			if len(my_current_threads) < self.threads:
				if thread_nr <= total_jobs:
					if thread_nr < total_jobs:
						new_thread = self.thread_list[thread_nr]
						new_thread.start()
						my_current_threads.append(new_thread)
					thread_nr += 1
					
			for thread in my_current_threads:
				if not thread.isAlive():
					print("finished thread: " + str(thread))
					self.write_result_line(thread)
					thread.set_handled()
			if len(my_current_threads) == 0:
				active_threads = False
			else:
				active_threads = True
			my_current_threads = [t for t in my_current_threads if not t.is_handled()]

	def calculate_DSS(self, max_anchor_boost, raw_dss_non_anchor, raw_dss_anchor, non_anchor_domains, anchor_domains):
		"""
		Calculates the domain similarity score based on the given information
		"""
		DSS_boosts = {}
		for anchor_boost in np.arange(1.0, float(1+max_anchor_boost), 1.0):
			if anchor_domains != 0 and non_anchor_domains != 0:
				# Calculate proportional weight to each kind of domain
				non_anchor_weight = float(non_anchor_domains) / (float(non_anchor_domains) + (float(anchor_domains) * float(anchor_boost)))
				anchor_weight = (float(anchor_domains) * float(anchor_boost)) / (float(non_anchor_domains) + (float(anchor_domains) * float(anchor_boost)))
				# scale the raw_dss_non_anchor and raw_dss_anchor
				DSS = (float(non_anchor_weight) * float(raw_dss_non_anchor)) + (float(anchor_weight) * float(raw_dss_anchor))
			elif anchor_domains == 0:
				DSS = raw_dss_non_anchor
			else: #only anchor domains present
				DSS = raw_dss_anchor
			DSS_boosts.update({anchor_boost:(DSS)}) # DSS is now still distance, which is needed (1 - DSS, to get similarity)
		return DSS_boosts

	def open_result_file(self):
		"""
		Opens the result file that will contain the final result
		"""
		self.rf = open(self.output_dir + "/connectivity_results.tsv", 'w')
		self.rf.write("File\tcutoff\tJw\tAIw\tDSSw\ttotal_bgcs\tnr_of_clusters\tbiggest_cluster_size\tsmallest_cluster_size\tnr_of_singletons\n")

	def write_result_line(self, finished_thread):
		"""
		Writes the results to the final output file
		"""
		print("writing results of: " + str(finished_thread))
		self.rf.write(finished_thread.get_result())


class NetworkConnectingThread(threading.Thread):
	"""
	A class to use for multithreading, extends the already existing threading.Thread
	"""
	def __init__(self, Jw, AIw, DSSw, distance_dict, current_file, ID):
		"""
		Initiates the class with the given arguments
		"""
		threading.Thread.__init__(self)
		self.ID = ID
		self.result_line = ""
		self.Jw = Jw
		self.AIw = AIw
		self.DSSw = DSSw
		self.distance_dict = distance_dict
		self.current_file = current_file
		self.handled = False

	def run(self):
		"""
		When Thread.start() is called on a Thread object, this will be called by the start() method in turn.
		Calculates the connectivity per file, with different cutoffs
		"""
		bgc_lists = dict()
		for bgc_pair in self.distance_dict.keys():
			bgc1 = bgc_pair[0]
			bgc2 = bgc_pair[1]
			Jaccard = self.distance_dict[bgc_pair]["J"]
			AI = self.distance_dict[bgc_pair]["AI"]
			DSS = self.distance_dict[bgc_pair]["DSS"]
			for anchor_boost in DSS.keys():
				if not bool(bgc_lists):
					bgc_lists = {key: {float(cutoff)/10:set() for cutoff in range(2, 9, 1)} for key in DSS.keys()}
				distance = (float(Jaccard) * self.Jw) + (float(AI) * self.AIw) + (float(DSS[anchor_boost]) * self.DSSw)
				for cutoff in range(2, 9, 1):
					cutoff = float(cutoff) / 10.0
					if distance < cutoff:
						bgc_lists[anchor_boost][cutoff].add((bgc1,bgc2))
					else:
						bgc_lists[anchor_boost][cutoff].add((bgc1, bgc1)) # only connected to itself
						bgc_lists[anchor_boost][cutoff].add((bgc2, bgc2)) # only connected to itself

		for ab in bgc_lists.keys():
			for cutoff in bgc_lists[ab].keys():
				self.current_bgcs = set(bgc_lists[ab][cutoff])
				self.connect_bgcs(cutoff)
				
		time_end = time.clock()

	def connect_bgcs(self, cutoff):
		"""
		Connects the bgcs in the list by merging sets that have a union
		"""
		g = networkx.Graph()
		for bgc_set in self.current_bgcs:
			g.add_edge(*bgc_set)
		self.create_result_line(networkx.connected_components(g) , cutoff)

	def create_result_line(self, result, cutoff):
		"""
		Writes the results to the final output file
		"""
		cluster_summaries = [len(x) for x in result]
		total_bgcs = sum(cluster_summaries)
		nr_of_singletons = cluster_summaries.count(1)
		nr_of_clusters = len(cluster_summaries) - nr_of_singletons
		try:
			biggest_cluster_size = max(cluster_summaries)
			if 1 in cluster_summaries:
				cluster_summaries.remove(1)
			smallest_cluster_size = min(cluster_summaries)
		except:
			biggest_cluster_size = 1
			smallest_cluster_size = 1
		self.result_line += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.current_file, cutoff, self.Jw, self.AIw, self.DSSw, total_bgcs, nr_of_clusters, biggest_cluster_size, smallest_cluster_size, nr_of_singletons)

	def get_result(self):
		"""
		Returns the results
		"""
		return self.result_line
	
	def set_handled(self):
		"""
		True if the thread is finished
		"""
		self.handled = True


	def is_handled(self):
		"""
		Returns if the thread is finished running
		"""
		return self.handled
		
def main():
	cc = ConnectivityCalculator()
	cc.getCmdArguments()
	cc.open_result_file()
	cc.parse_network_files()
	cc.prepare_distance_with_weights_calculations()

if __name__ == "__main__":
	main()
