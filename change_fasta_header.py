import sys, os

def main():
	new_file = open("../Genomes/Streptomyces_dna.faa", 'w')
	for infile in sys.argv[1:]:
		for line in open(infile):
			if line.startswith(">"):
				new_file.write(">" + line.strip().split("|")[1].split(".")[0] + line.strip().split("_")[2].split(".")[0] + "\n")
			else:
				new_file.write(line)
	
	
if __name__ == "__main__":
	main()
