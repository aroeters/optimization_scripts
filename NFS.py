#! /usr/bin/python3
######################################
# Author: Arne Roeters				 #
# MSc Thesis						 #
# Bioinformatics group wageningen	 #
# Purpose:							 #
# Scoring network files on a per GCF #
# basis where the sum of all scores  #
# is the total score for the file    #
######################################

## Imports
import argparse
import sys
import os
import bz2
from math import log2

class NetworkScoreCalculator:
	def __init__(self, network_file, AP_file, bgc_to_compound_file, group_file, output_dir):
		'''
		Initializer of the class NetworkScoreCalculator
		'''
		if output_dir.endswith("/"):
			self.out = output_dir
		else:
			self.out = output_dir + "/"
		self.nf = network_file
		self.apf = AP_file
		self.btc = bgc_to_compound_file
		self.gf = group_file
		self.AP_groups = dict()
		self.gcf_links = dict()
		self.gcf_scores = dict()
		self.total_score = float(0.00)
		self.total_nr_bgcs = 0
		self.bgc_to_compound = dict()
		self.compound_to_group = dict()
		self.groups_per_gcf = dict()
		self.edges_per_gcf = dict()
	
	def get_compound_to_class(self):
		'''
		Parses the group file that has the groups corresponding to the compounds.
		'''
		with open(self.gf) as f:
			next(f)
			for line in f:
				compound_class, compound = [line.strip().split("\t")[i] for i in [2,3]]
				self.compound_to_group.update({compound:compound_class})
	
	def get_bgc_to_compounds(self):
		'''
		Parses the bgc to compounds file to create a dict that translates the BGC to the corresponding compound
		'''
		for line in open(self.btc):
			bgc, compound = line.strip().split("\t")
			if compound in self.compound_to_group.keys():
				if bgc in self.bgc_to_compound:
					# To check if the same bgc does not have two classifications
					if self.compound_to_group[compound] != self.compound_to_group[self.bgc_to_compound[bgc]]:
						print("Double bgc entry with different classification for:\n" + bgc + ":" + compound + ":" + self.compound_to_group[compound]+"\nWith: " + self.bgc_to_compound[bgc] + ":" + self.compound_to_group[self.bgc_to_compound[bgc]]+"\n")
				else:
					self.bgc_to_compound.update({bgc:compound})
			
	def get_affinity_propagation_groups(self):
		'''
		Reads in the affinity propagation groups and adds them to the dict for later processing
		'''
		try:
			for line in open(self.apf):
				split_line = line.strip().split(",")
				self.AP_groups.update({split_line[0]:split_line[1::]})
				self.groups_per_gcf.update({split_line[0]:set()})
				self.gcf_scores.update({split_line[0]:0})
		except IOError as e:
			print("Failed to open propagation file: " + self.apf + "\nError:\n")
			sys.exit(e) 

	def get_all_clusters_from_network(self):
		'''
		Retrieves the information from the network file
		'''
		try:
			if self.nf.endswith(".bz2"):
				network_file = bz2.BZ2File(self.nf, 'r')
			else:
				network_file = open(self.nf)
			with network_file as f:
				starting_column = f.readline().split("\t").index("Raw distance")
				for line in f:
					bgc1, bgc2, distance = [line.split("\t")[i] for i in [0, 1, starting_column]]
					try:
						compound = self.bgc_to_compound[bgc1]
						group1 = self.compound_to_group[compound]
					except KeyError as e:
						group1 = "unknown"
					try:
						compound = self.bgc_to_compound[bgc2]
						group2 = self.compound_to_group[compound]
					except KeyError as e:
						group2 = "unknown"
					
					if bgc1 not in self.gcf_links.keys():
						self.gcf_links.update({bgc1:{bgc2:[float(float(1.00)-float(distance)), group1, group2]}})
					else:
						self.gcf_links[bgc1].update({bgc2:[float(float(1.00)-float(distance)), group1, group2]})
		except IOError as e:
			print("Failed to open network file: " + self.nf + "\nError:\n")
			sys.exit(e)
	
	def calculate_cluster_score(self):
		'''
		Calculates the score of every gcf
		'''
		for gcf in self.AP_groups.keys():
			self.edges_per_gcf.update({gcf:{"b":0, "w":0, "n":0}}) # b = between edge, w = within edge, n = neutral edge
			# First calculate how many different groups are in the GCF
			for bgc in self.AP_groups[gcf]:
				self.total_nr_bgcs += 1
				try:
					group = self.compound_to_group[self.bgc_to_compound[bgc]]
					self.groups_per_gcf[gcf].add(group)
				except KeyError as e:
					continue
			
			# for every combination of bgc's in the gcf where bgc1 != bgc2
			for bgc1 in self.AP_groups[gcf]:
				for bgc2 in self.AP_groups[gcf]:
					if bgc1 != bgc2:
						try:
							# if the link exists between bgc 1 and 2
							similarity, group1, group2 = self.gcf_links[bgc1][bgc2]
							if group1 == "unknown" or group2 == "unknown":
								self.gcf_scores[gcf] += float(similarity)
								self.edges_per_gcf[gcf]["n"] += 1
							elif group1 != group2:
								#~ self.gcf_scores[gcf] -= ((float(1.00) + float(similarity))*len(self.groups_per_gcf[gcf]))
								self.gcf_scores[gcf] -= (float(1.00) + float(similarity)) # try with lower penalty
								self.edges_per_gcf[gcf]["b"] += 1
							else:
								self.gcf_scores[gcf] += (float(1.00) + float(similarity))
								self.edges_per_gcf[gcf]["w"] += 1
						except KeyError as e:
							continue
			self.total_score += self.gcf_scores[gcf]
	
	def calculate_cluster_score_entropy(self):
		'''
		Calculates the score of every gcf using the entropy of that given gcf
		'''
		for gcf in self.AP_groups.keys():
			entropy_dict = {}
			total_bgcs_of_gcf = len(self.AP_groups[gcf])
			self.edges_per_gcf.update({gcf:{"b":0, "w":0, "n":0, 't':0}}) # b = between edge, w = within edge, n = neutral edge
			
			for bgc in self.AP_groups[gcf]:
				self.total_nr_bgcs += 1
				try:
					group = self.compound_to_group[self.bgc_to_compound[bgc]]
					if not group in entropy_dict:
						entropy_dict.update({group:1})
					else:
						entropy_dict[group] += 1
				except KeyError as e:
					if "unknown" in entropy_dict:
						entropy_dict["unknown"] += 1
					else:
						entropy_dict.update({"unknown":1})
			# gets the entropy of the gcf
			entropy_of_gcf = sum([-1 * (float(entropy_dict[x]/total_bgcs_of_gcf)*log2(float(entropy_dict[x]/total_bgcs_of_gcf))) for x in entropy_dict.keys()])
			
			# for every combination of bgc's in the gcf where bgc1 != bgc2
			for bgc1 in self.AP_groups[gcf]:
				for bgc2 in self.AP_groups[gcf]:
					if bgc1 != bgc2:
						try:
							# if the link exists between bgc 1 and 2
							similarity, group1, group2 = self.gcf_links[bgc1][bgc2]
							if group1 == "unknown" or group2 == "unknown":
								self.gcf_scores[gcf] += float(similarity)
								self.edges_per_gcf[gcf]["n"] += 1
							elif group1 != group2:
								self.gcf_scores[gcf] -= (float(1.00) + float(similarity)) * entropy_of_gcf
								self.edges_per_gcf[gcf]["b"] += 1
							else:
								self.gcf_scores[gcf] += (float(1.00) + float(similarity))
								self.edges_per_gcf[gcf]["w"] += 1
						except KeyError as e:
							continue
			self.total_score += self.gcf_scores[gcf]
	
	def write_results_to_file(self):
		'''
		Writes the final scores to the result file
		'''
		# to create more consistency in the file names
		nf_name = self.apf.replace("\\", "/").split("/")[-1].replace(".csv", "").replace("affinity_propagation_","").replace("NRP_Hybrids", "NRP-Hybrids").replace("all_mix", "all-mix").replace("all_NRPS", "allNRPS")
		
		result_file = open("{}{}_score_entropy.tsv".format(self.out, nf_name), "w")
		result_file.write("Total_score_normalized:{}\nTotal_score:{}\nTotal_nr_of_bgcs:{}\n\nGcf\tNr_of_bgcs\tScore\tScore_normalized_on_edges\tWithin\tNeutral\tBetween\n".format(str(self.total_score/self.total_nr_bgcs), str(self.total_score), str(self.total_nr_bgcs)))
		
		for gcf in self.gcf_scores.keys():
			if self.edges_per_gcf[gcf]["w"] + self.edges_per_gcf[gcf]["b"] + self.edges_per_gcf[gcf]["n"] != 0:
				normalized_score = str(float(self.gcf_scores[gcf])/float(self.edges_per_gcf[gcf]["w"] + self.edges_per_gcf[gcf]["b"] + self.edges_per_gcf[gcf]["n"]))
			else:
				normalized_score = "0"
			result_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gcf, str(len(self.AP_groups[gcf])), str(self.gcf_scores[gcf]), normalized_score, self.edges_per_gcf[gcf]["w"], self.edges_per_gcf[gcf]["n"], self.edges_per_gcf[gcf]["b"]))
	
	def create_filtered_network(self):
		'''
		Creates a filtered network file based on the affinity propagation
		'''
		filtered_network_file = open(self.out + self.apf.replace("\\", "/").split("/")[-1].replace(".csv", "").replace(".bz2", "").replace("affinity_propagation_","") + "_AP_filtered.network", 'w')
		possible_combinations = dict()
		#Get all possible AP combinations possible
		for gcf in self.AP_groups.keys():
			for bgc1 in self.AP_groups[gcf]:
				for bgc2 in self.AP_groups[gcf]:
					if bgc1 != bgc2:
						if bgc1 not in possible_combinations.keys():
							possible_combinations.update({bgc1:{bgc2:gcf}})
						else:
							possible_combinations[bgc1].update({bgc2:gcf})
						if bgc2 not in possible_combinations.keys():
							possible_combinations.update({bgc2:{bgc1:gcf}})
						else:
							possible_combinations[bgc2].update({bgc1:gcf})
		
		with open(self.nf) as f:
			filtered_network_file.write("GCF\t" + f.readline())
			for line in f:
				bgc1, bgc2 = [line.split("\t")[i] for i in [0,1]]
				if any([x in possible_combinations.keys() for x in [bgc1, bgc2]]):
					try:
						if bgc2 in possible_combinations[bgc1]:
							filtered_network_file.write(possible_combinations[bgc1][bgc2] + "\t" + line)
					except KeyError:
						if bgc1 in possible_combinations[bgc2]:
							filtered_network_file.write(possible_combinations[bgc2][bgc1] + "\t" + line)
							
			
def get_cmd_arguments():
	parser = argparse.ArgumentParser(prog='NFS.py')
	parser.add_argument("-a", "--affinity_propagation_file", help="The affinity propagation file\n(obtained from the script Affinity_Propagation_new_weights.py)\nwith complete path to the file", required=True)
	parser.add_argument("-b", "--compound_file", help="The bgc to compound file, with complete path to the file",required=True)
	parser.add_argument("-g", "--group_file", help="The compound to class file, with complete path to the file",required=True)
	parser.add_argument("-n", "--network_file", help="The BGC network file, with complete path to the file",required=True)
	parser.add_argument("-o", "--output_dir", help="The directory to put the resulting score file in.", required=False, default="./")
	parser.add_argument("-f", "--filtered_network", help="If used, the program will also output a filtered network file based on the affinity propagation.", action='store_true')
	args = parser.parse_args()
	
	if not os.path.isfile(args.network_file) or not os.path.isfile(args.affinity_propagation_file) or not os.path.isfile(args.compound_file):
		sys.exit("Check your input files for the following:\n\t- Correct full path given\n\t- Correct name\nGiven:\n\t{}\n\t{}\n\t{}".format(args.network_file, args.affinity_propagation_file, args.compound_file))
	elif not os.path.isdir(args.output_dir):
		sys.exit("The specified directory:\n\n{}\n\ndoes not exist.".format(args.output_dir))
	else:
		return args.network_file, args.affinity_propagation_file, args.compound_file, args.group_file, args.output_dir, args.filtered_network

def main():
	network_file, AP_file, compound_file, group_file, output_dir, filtered_network = get_cmd_arguments()
	NSC = NetworkScoreCalculator(network_file, AP_file, compound_file, group_file, output_dir)
	NSC.get_compound_to_class()
	NSC.get_bgc_to_compounds()
	NSC.get_all_clusters_from_network()
	NSC.get_affinity_propagation_groups()
	#~ NSC.calculate_cluster_score()
	NSC.calculate_cluster_score_entropy()
	NSC.write_results_to_file()
	if filtered_network:
		NSC.create_filtered_network()

if __name__ == "__main__":
	main()
