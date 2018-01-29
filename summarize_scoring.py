import sys
import glob

def summarize_scoring(out, last_part_of_file):
	"""
	Goes through all scoring files that are in the given output folder 
	and summarizes a part of the output
	""" 
	file_path = out+"/results_scoring.txt"
	final_results = open(file_path, 'w')
	print("Final results in: " + file_path)
	final_results.write("file\tscore\tBGCs\tnormalized_score\tBGCs_max_GCF\tmax_GCF_score\tmax_GCF_nr\tBGCs_min_GCF\tmin_GCF_score\tmin_GCF_nr\tJI\tAI\tDSS\tAB\tcutoff\n")
	scoring_files = glob.glob(out + "/*"+last_part_of_file)
	print(len(scoring_files))
	
	
	max_GCF_score = float('-inf')
	max_GCF = ""
	nr_of_BGCs_max = 0
	min_GCF_score = float('inf')
	min_GCF = ""
	nr_of_BGCs_min = 0
	

	for sf in sorted(scoring_files):
		total_edge_based_score = 0.0
		total_edge_lines = 0.0
		final_results.write(sf.split("/")[-1].replace(last_part_of_file, "\t"))
		cutoff, JI, AI, DSS, AB = sf.split("/")[-1].replace("new_", "").replace("_score.txt", "").split("_")[1:6]
		total_score_info = True
		for line in open(sf):
			if total_score_info:
				if line.startswith("Gcf"):
					total_score_info = False
				else:
					if line != "\n" and "normalized" not in line:
						info = round(float(line.split(":")[-1].strip()), 3)
						final_results.write(str(info)+"\t")
			else:
				gcf, nr_of_BGC, score, edge_based_score = [line.split("\t")[i].strip() for i in [0,1,2,3]]
				if float(score) > max_GCF_score:
					max_GCF_score = float(score)
					max_GCF = gcf
					nr_of_BGCs_max = nr_of_BGC
				elif float(score) < min_GCF_score:
					min_GCF_score = float(score)
					min_GCF = gcf
					nr_of_BGCs_min = nr_of_BGC
				total_edge_based_score += float(edge_based_score)
				total_edge_lines += 1
		
		if total_edge_lines == 0:
			normalized = 0.0
			nr_of_BGCs_max = 0 
			gcf_max = "1.0"
			max_GCF = "0"
			nr_of_BGCs_min = 0
			gcf_min = "1.0"
			min_GCF = "0"
			final_results.write("-500.0\t-500.0\t")
			
		else:
			normalized = round(total_edge_based_score/total_edge_lines,3)
			gcf_max = str(round(max_GCF_score, 3))
			gcf_min = str(round(min_GCF_score, 3))
			
		final_results.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(normalized ,nr_of_BGCs_max, gcf_max, max_GCF, nr_of_BGCs_min, gcf_min, min_GCF, JI.replace("Jw", ""), AI.replace("AIw",""), DSS.replace("DSSw", ""), AB.replace("AB",""), cutoff.replace("c", "")))
	final_results.close()

def main():
	outdir = sys.argv[1]
	last_part_of_file = sys.argv[2]
	summarize_scoring(outdir, last_part_of_file)

if __name__ == "__main__":
	main()
