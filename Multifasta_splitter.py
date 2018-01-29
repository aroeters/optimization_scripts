################################
##  Author: Arne Roeters      ##
##  Splits a single multi     ##
##  fasta file from the ncbi  ##
##  database into multiple    ##
##  fasta files per organism  ##
################################ 

import os, sys, argparse

def getCMDArguments():
	"""
	Retrieves the arguments from the cmd
	"""
	parser = argparse.ArgumentParser(prog='Multifasta_splitter.py')
	
	parser.add_argument("-f" ,"--faa", help="the mutlifasta file to split", required=True)
	parser.add_argument("-o" ,"--out", help="Output directory for the fasta files", required=True)
	parser.add_argument("-s" ,"--shotgun_data", help="Flag if the data comes from whole genome shotgun sequencing", required=False, action='store_true')
	parser.add_argument("-n" ,"--names", help="A file with accession and organism separated by a tab", required=False, default='/home/roete009/Thesis/python/temp/Streptomyces_names.tsv')
	
	args = parser.parse_args()
	
	faa = args.faa
	out = args.out
	shotgun = args.shotgun_data
	
	if not os.path.isfile(faa):
		sys.exit("Please specify a multifasta file\nGiven : {}".format(faa))
	if not os.path.isdir(out):
		sys.exit("Please specify a valid output directory\nGiven : {}".format(out))
	if not out.endswith("/"):
		out = out + "/"
	return faa, out, shotgun, args.names
	
def splitFaaFile(faa, out, names):
	"""
	Splits the given faa file into multiple single organism files
	"""
	current_organism = "ThisIsNotAnID"		
	accessionToOrganism = getSpeciesFromAccession(names)
	for line in open(faa):
		if not line.strip() == "":
			if line.startswith(">"):
				if current_organism not in line and not "NZ_" in line:
					current_organism = line.split()[0].split("|")[1].split(".")[0]
					#creates a fasta file with the organism id and .fa as extension
					output_file = open(out + current_organism.strip() + ".faa" , 'w')
					# Writes the first header to the new file
					output_file.write(line)
				else:
					output_file.write(line)
			else:
				output_file.write(line)

def splitShotgunToFaa(faa, out, names):
	"""
	Splits the contigs in the shotgun data into files with a single organism
	"""
	accessionToOrganism = getSpeciesFromAccession(names)
	organism_start = "initial_start"
	for line in open(faa):
		if not line.strip() == "":
			if line.startswith(">"):
				if organism_start not in line:
					organism_start = line.split()[0].split("|")[1].split(".")[0][0:8]
					#creates a fasta file with the organism id and .fa as extension
					output_file = open(out + organism_start + ".faa" , 'w')
					# Writes the first header to the new file
					output_file.write(line)
				else:
					output_file.write(line)
			else:
				output_file.write(line)
	
def getSpeciesFromAccession(infile):
	"""
	Creates a dictionary that converts accession to organism name
	"""
	accessionToOrganism = dict()
	for line in open(infile):
		accession, organism = line.split("\t")
		accessionToOrganism.update({accession:organism})
	return accessionToOrganism
			
def main():
	faa, out, shotgun, names = getCMDArguments()
	if not shotgun:
		splitFaaFile(faa, out, names)
	else:
		splitShotgunToFaa(faa, out, names)
	
if __name__ == "__main__":
	main()
