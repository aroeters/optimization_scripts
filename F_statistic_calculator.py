import glob
import os
import sys
import argparse

class FStatisticCalculator():
    def __init__(self):
        """
        Initiates the class with the created variables
        """
        self.gf = "" # group file
        self.btc = "" # BGC to compound file
        self.output_dir = ""
        self.ap_files = list() # affinity propagation files
        self.bgc_to_compound = dict()
        self.compound_to_group = dict()
        self.group_size = dict()
        self.dominant_group_counts = dict()
        self.complete_network_score = 0
        
    def getCmdArguments(self):
        """
        Retrieves the command line arguments
        """
        parser = argparse.ArgumentParser(prog='F_statistic_calculator.py')
        parser.add_argument("-a", "--affinity_propagation_dir", help="The directory with the to be used affinity propagation files in it", required=True)
        parser.add_argument("-b", "--compound_file", help="The bgc to compound file, with complete path to the file",required=False, default="/home/roete009/Thesis/python/temp/bgcs_to_compounds_ALL.tsv")
        parser.add_argument("-g", "--group_file", help="The compound to class file, with complete path to the file",required=False, default="/home/roete009/Thesis/python/temp/training_set_rev03_classes.tsv")
        parser.add_argument("-o", "--output_dir", help="The directory to output the final results in", required=False, default="./")
        args = parser.parse_args()
        
        if os.path.isfile(args.group_file):
            self.gf = args.group_file
        else:
            sys.exit("The given group file is not found")
        
        if os.path.isfile(args.compound_file):
            self.btc = args.compound_file
        else:
            sys.exit("The given compound file is not found")
        
        if os.path.isdir(args.output_dir):
            self.output_dir = args.output_dir
        else:
            sys.exit("The given output directory is not valid")
            
        if os.path.isdir(args.affinity_propagation_dir):
            files = glob.glob(args.affinity_propagation_dir + "/*affinity_propagation*")
            if len(files) != 0:
                self.ap_files = files
            else:
                sys.exit("No affinity propagation files found in the given directory")
        else:
            sys.exit("The given AP directory is not valid")
        print("Command line arguments retrieved...")
        
    def get_compound_to_group(self):
        '''
        Parses the group file that has the groups corresponding to the compounds.
        '''
        with open(self.gf) as f:
            next(f)
            for line in f:
                compound_group, compound = [line.strip().split("\t")[i] for i in [2,3]]
                self.compound_to_group.update({compound:compound_group})
                if compound_group not in self.group_size.keys():
                    self.group_size.update({compound_group:0})
                self.group_size[compound_group] += 1
        print("Compound to group processed...")
    
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
        print("BGC to compound processed...")
    
    def convert_to_groups(self, gcf_bgcs):
        """
        Converts the BGCS to groups keeping the unknown groups out
        """
        gcf_groups = list()
        for bgc in gcf_bgcs:
            try:
                gcf_groups.append(self.compound_to_group[self.bgc_to_compound[bgc]])
            except: # if not found in either of the lists
                continue # do nothing when the bgc group is unknown    
        return gcf_groups                    
    
    def parse_ap_files(self):
        """
        parses the affinity propagation files
        """
        print("Processing all affinity propagation files")
        for ap_file in self.ap_files:
        # create a dictionary with all group names as keys and an empty list as the value
            self.dominant_group_counts = {group:[0, 0] for group in self.compound_to_group.values()}
            self.complete_network_score = 0
            for line in open(ap_file):
                splitline = line.split(",")
                gcf_nr = splitline[0]
                gcf_bgcs = splitline[1:]
                gcf_groups = self.convert_to_groups(gcf_bgcs) # Also removes the unknown groups
                self.process_group_dominance(self.get_group_counts(gcf_groups))
            self.calculate_network_score()
            self.write_result_to_file(ap_file)
        
    def process_group_dominance(self, group_counts):
        """
        Retrieves the gcf of all GCFs
        """
        total = sum(group_counts.values())
        for group in group_counts.keys():
            if self.dominant_group_counts[group][0] < group_counts[group]:
                self.dominant_group_counts[group] = [group_counts[group], total]
            elif self.dominant_group_counts[group][0] == group_counts[group]:
                # This is done to penalize for groups being equally divided over multiple GCFs
                if total > self.dominant_group_counts[group][1]:
                    self.dominant_group_counts[group][1] = total
            
    def get_group_counts(self, gcf_groups):
        """
        For every compound group the number of occurences in the gcf is counted
        """
        group_counts = dict()
        # count all occurences for each individual group that is found within the gcf
        for group in set(gcf_groups):
            group_counts.update({group:gcf_groups.count(group)})
        return group_counts    
        
    def calculate_network_score(self):
        """
        Calculates the F-statistic based on the clusters in the GCF based on the dominant group in the GCF
        """
        total_cluster = 0
        for group in self.dominant_group_counts:
            total_cluster += 1
            group_count, total = self.dominant_group_counts[group]

            if total != 0:
                # The ratio of the total number of nodes in the curated group vs the number of nodes
                # in the GCF with the largest number of members in that group
                completeness = float(group_count) / float(self.group_size[group])
                # The purity determines how many other BGCs that do not belong to the dominant group are in the GCF
                purity = float(group_count) / float(total)
                # *0.5 to get a score between 0.0 and 1.0
                self.complete_network_score += (completeness * 0.50) + (purity * 0.50)
        self.complete_network_score = self.complete_network_score / total_cluster
    
    def open_output_file(self):
        """
        Creates the output file
        """
        self.file_out = open(self.output_dir + "/f_statistic_scores.tsv", 'w')
        
    def write_result_to_file(self, ap_file):
        """
        Writes the results for all the AP files to a final results file in a tsv format
        """
        ap_file = ap_file.split("affinity_propagation_")[1]
        cutoff = ap_file.split("_c")[1].split("_")[0]
        jaccard = ap_file.split("_Jw")[1].split("_")[0]
        adjacency = ap_file.split("_AIw")[1].split("_")[0]
        DSS = ap_file.split("_DSSw")[1].split("_")[0]
        AB = ap_file.split("_AB")[1].split(".csv")[0]
        self.file_out.write(ap_file + "\t" + str(self.complete_network_score) + "\t" + cutoff + "\t" + jaccard + "\t" + adjacency + "\t" + DSS + "\t" + AB + "\n")
        
def main():
    F = FStatisticCalculator()
    F.getCmdArguments()
    F.get_compound_to_group()
    F.get_bgc_to_compounds()
    F.open_output_file()
    F.parse_ap_files()
    
if __name__ == "__main__":
    main()
