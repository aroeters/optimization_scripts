from scipy.sparse import lil_matrix
import pysapc
import numpy as np
import bz2
import argparse
import os, sys
from sklearn.cluster import AffinityPropagation
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import pam250 as scoring_matrix



global bgc_fasta_folder
global pfd_folder
##########################################################
## Needed to create the json file for the visualization ##
##########################################################
#~ def get_gbk_files(inputdir, bgc_fasta_folder, min_bgc_size, exclude_gbk_str, bgc_info):
def get_gbk_files(inputdir, min_bgc_size, bgc_info):
    """Searches given directory for genbank files recursively, will assume that
    the genbank files that have the same name are the same genbank file. 
    Returns a dictionary that contains the names of the clusters found as keys
    and a list that contains [0] a path to the genbank file and [1] the 
    samples that the genbank file is a part of.
    Extract and write the sequences as fasta files if not already in the Fasta 
    folder.
    return: {cluster_name:[genbank_path,[s_a,s_b...]]}
    """

    genbankDict = {}
    file_counter = 0
    processed_sequences = 0
    biosynthetic_genes = set()
    product_list_per_record = []
    fasta_data = []
    save_fasta = True
    contig_edge = False
    
    current_dir = ""
    for dirpath, dirnames, filenames in os.walk(inputdir):
        head, tail = os.path.split(dirpath)

        if current_dir != tail:
            current_dir = tail

        genbankfilelist = []

        for fname in filenames:
            if fname[-3:] != "gbk":
                continue
            
            clusterName = fname[:-4]
            
            if "_ORF" in fname:
                print(" Skipping file " + fname + " (string '_ORF' is used internally)")
                continue
            
            if " " in fname:
                sys.exit("\nError: Input GenBank files should not have spaces in their filenames as HMMscan cannot process them properly ('too many arguments').")
            
            # See if we need to write down the sequence
            outputfile = os.path.join(bgc_fasta_folder, clusterName + '.fasta')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                save_fasta = False
            else:
                save_fasta = True
                
            try:
                # basic file verification. Substitutes check_data_integrity
                records = list(SeqIO.parse(os.path.join(dirpath,fname), "genbank"))
            except ValueError as e:
                print("   Error with file " + os.path.join(dirpath, fname) + ": \n    '" + str(e) + "'")
                print("    (This file will be excluded from the analysis)")
                continue
            else:
                bgc_size = 0
                cds_ctr = 0
                product = "no type"
                del product_list_per_record[:]
                
                max_width = 0 # This will be used for the SVG figure
                record_count = 0
                
                for record in records:
                    record_count += 1
                    bgc_size += len(record.seq)
                    if len(record.seq) > max_width:
                        max_width = len(record.seq)
                    
                    for feature in record.features:
                        if "cluster" in feature.type and "product" in feature.qualifiers:
                            if len(feature.qualifiers["product"]) > 1:
                                print("  WARNING: more than product annotated in record " + str(record_count) + ", " + fname)
                                break
                            else:
                                product_list_per_record.append(feature.qualifiers["product"][0].replace(" ",""))
                
                        # Get biosynthetic genes + sequences
                        if feature.type == "CDS":
                            cds_ctr += 1
                                    
                            CDS = feature
                            gene_id = ""
                            if "gene" in CDS.qualifiers:
                                # In principle, we should keep a list of genes with
                                # the same id (isoforms) and only keep the largest
                                # TODO
                                gene_id = CDS.qualifiers.get('gene',"")[0]
                            
                            protein_id = ""
                            if "protein_id" in CDS.qualifiers:
                                protein_id = CDS.qualifiers.get('protein_id',"")[0]
                            
                            # nofuzzy_start/nofuzzy_end are obsolete
                            # http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html#nofuzzy_start
                            gene_start = max(0, int(CDS.location.start))
                            gene_end = max(0, int(CDS.location.end))
                            direction = CDS.location.strand
                            
                            if direction == 1:
                                strand = '+'
                            else:
                                strand = '-'
                                
                            fasta_header = clusterName + "_ORF" + str(cds_ctr)+ ":gid:" + str(gene_id) + ":pid:" + str(protein_id) + ":loc:" + str(gene_start) + ":" + str(gene_end) + ":strand:" + strand
                            fasta_header = fasta_header.replace(">","") #the coordinates might contain larger than signs, tools upstream don't like this
                            fasta_header = fasta_header.replace(" ", "") #the domtable output format (hmmscan) uses spaces as a delimiter, so these cannot be present in the fasta header

                            if "sec_met" in feature.qualifiers:
                                if "Kind: biosynthetic" in feature.qualifiers["sec_met"]:
                                    biosynthetic_genes.add(fasta_header)

                            fasta_header = ">"+fasta_header
                            

                            if save_fasta:
                                if 'translation' in CDS.qualifiers.keys():
                                    prot_seq = CDS.qualifiers['translation'][0]
                                # If translation isn't available translate manually, this will take longer
                                else:
                                    nt_seq = CDS.location.extract(record.seq)
                                    
                                    # If we know sequence is an ORF (like all CDSs), codon table can be
                                    #  used to correctly translate alternative start codons.
                                    #  see http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc25
                                    # If the sequence has a fuzzy start/end, it might not be complete,
                                    # (therefore it might not be the true start codon)
                                    # However, in this case, if 'translation' not available, assume 
                                    #  this is just a random sequence 
                                    complete_cds = False 
                                    
                                    # More about fuzzy positions
                                    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc39
                                    fuzzy_start = False 
                                    if str(CDS.location.start)[0] in "<>":
                                        complete_cds = False
                                        fuzzy_start = True
                                        
                                    fuzzy_end = False
                                    if str(CDS.location.end)[0] in "<>":
                                        fuzzy_end = True
                                    
                                    #for protein sequence if it is at the start of the entry assume 
                                    # that end of sequence is in frame and trim from the beginning
                                    #if it is at the end of the genbank entry assume that the start 
                                    # of the sequence is in frame
                                    reminder = len(nt_seq)%3
                                    if reminder > 0:
                                        if fuzzy_start and fuzzy_end:
                                            print("Warning, CDS (" + clusterName + ", " + CDS.qualifiers.get('locus_tag',"")[0] + ") has fuzzy start and end positions, and a sequence length not multiple of three. Skipping")
                                            break
                                        
                                        if fuzzy_start:
                                            if reminder == 1:
                                                nt_seq = nt_seq[1:]
                                            else:
                                                nt_seq = nt_seq[2:]
                                        # fuzzy end
                                        else:
                                            #same logic reverse direction
                                            if reminder == 1:
                                                nt_seq = nt_seq[:-1]
                                            else:
                                                nt_seq = nt_seq[:-2]
                                    
                                    # The Genetic Codes: www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
                                    if "transl_table" in CDS.qualifiers.keys():
                                        CDStable = CDS.qualifiers.get("transl_table", "")[0]
                                        prot_seq = str(nt_seq.translate(table=CDStable, to_stop=True, cds=complete_cds))
                                    else:
                                        prot_seq = str(nt_seq.translate(to_stop=True, cds=complete_cds))
                                        
                                fasta_data.append((fasta_header, prot_seq))
                
                if bgc_size > min_bgc_size:  # exclude the bgc if it's too small
                    file_counter += 1
                    # check what we have product-wise
                    # In particular, handle different products for multi-record files
                    product_set = set(product_list_per_record)
                    if len(product_set) == 1: # only one type of product
                        product = product_list_per_record[0]
                    elif "other" in product_set: # more than one, and it contains "other"
                        if len(product_set) == 2:
                            product = list(product_set - set(['other']))[0] # product = not "other"
                        else:
                            product = "-".join(product_set - set(['other'])) # likely a hybrid
                    else:
                        product = "-".join(product_set) # likely a hybrid
                     
                    
                    # assuming that the definition field is the same in all records
                    # product: antiSMASH predicted class of metabolite
                    # gbk definition
                    # number of records (for Arrower figures)
                    # max_width: width of the largest record (for Arrower figures)
                    # id: the GenBank's accession
                    #bgc_info[clusterName] = (product, records[0].description, len(records), max_width, records[0].id, biosynthetic_genes.copy())
                    # TODO contig_edge annotation is not present for antiSMASH v < 4
                    # Perhaps we can try to infer if it's in a contig edge: if
                    # - first biosynthetic gene start < 10kb or
                    # - max_width - last biosynthetic gene end < 10kb (but this will work only for the largest record)
                    bgc_info[clusterName] = bgc_data(records[0].id, records[0].description, product, len(records), max_width, biosynthetic_genes.copy(), contig_edge)

                    # TODO why re-process everything if it was already in the list?
                    # if name already in genbankDict.keys -> add current_dir
                    # else: extract all info
                    if clusterName in genbankDict.keys():
                        # Name was already in use. Use current_dir as the new sample's name
                        genbankDict[clusterName][1].add(current_dir) 
                    else:
                        # location of first instance of the file is genbankDict[clustername][0]
                        genbankDict.setdefault(clusterName, [os.path.join(dirpath, fname), set([current_dir])])
                        
                        # See if we need to write down the sequence
                        outputfile = os.path.join(bgc_fasta_folder, clusterName + '.fasta')
                        if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                            processed_sequences += 1
                        else:
                            with open(outputfile,'w') as fastaHandle:
                                for header_sequence in fasta_data:
                                    fastaHandle.write('%s\n' % header_sequence[0])
                                    fastaHandle.write('%s\n' % header_sequence[1])             
                else:
                    print(" Discarding " + clusterName +  " (size less than " + str(min_bgc_size) + " bp, was " + str(bgc_size) + ")")
                
                del fasta_data[:]
                biosynthetic_genes.clear()

    if file_counter == 0:
        sys.exit("\nError: There are no files to process")
    if file_counter == 1:
        sys.exit("\nError: Only one file found. Please input at least two files")

    print("\n Starting with " + str(file_counter) + " files")
    print(" Files that had its sequence extracted: " + str(file_counter - processed_sequences))

    return genbankDict

class bgc_data:
	def __init__(self, accession_id, description, product, records, max_width, biosynthetic_genes, contig_edge):
		# These two properties come from the genbank file:
		self.accession_id = accession_id
		self.description = description
		# AntiSMASH predicted class of compound:
		self.product = product
		# number of records in the genbank file (think of multi-locus BGCs):
		self.records = records
		# length of largest record (it will be used for ArrowerSVG):
		self.max_width = int(max_width)
		# Internal set of tags corresponding to genes that AntiSMASH marked 
		# as "Kind: Biosynthetic". It is formed as
		# clusterName + "_ORF" + cds_number + ":gid:" + gene_id + ":pid:" + protein_id + ":loc:" + gene_start + ":" + gene_end + ":strand:" + {+,-}
		self.biosynthetic_genes = biosynthetic_genes
		# AntiSMASH 4+ marks BGCs that sit on the edge of a contig
		self.contig_edge = contig_edge


class NetworkFileReader():
	def __init__(self, ibgc_file, Jw, AIw, DSSw, anchor_boost, output_dir, cutoff, used_ap, gbk_file_dir):
		"""
		Initializer of the class NetworkFileReader
		"""
		self.ibgc = ibgc_file
		self.DSSw = DSSw
		self.AIw = AIw
		self.Jw = Jw
		self.anchor_boost = anchor_boost
		self.bgcs = set()
		self.similarity_dict = {}
		if output_dir.endswith("/"):
			self.output_dir = output_dir
		else:
			self.output_dir = output_dir + "/"
		self.cutoff = float(cutoff)
		self.nf_name = self.ibgc.replace("\\", "/").split("/")[-1].replace(".network", "").replace(".bz2", "")
		self.used_ap = used_ap
		self.gbk_file_dir = gbk_file_dir

	def process_network_file(self):
		"""
		processes the network file line by line
		"""
		header=True
		starting_column=0
		if "bz2" in self.ibgc:
			networkFile = bz2.BZ2File(self.ibgc, 'r')
		else:
			networkFile = open(self.ibgc, 'r')
		for line in networkFile:
			if header:
				# This is done because there are two versions of the file, the one with only MIBIG misses some columns in between
				# It can now process both files
				starting_column = line.split("\t").index("Raw distance")  #skip header
				header=False
			else:
				splitLine = line.split("\t")
				# Jaccard, DSS, AI, nr_non_anchor_domains, nr_anchor_domains, DDS_non_anchor, DDS_anchor
				distance = self.process_network_line(splitLine[starting_column+2],splitLine[starting_column+3],splitLine[starting_column+4], splitLine[starting_column+7], splitLine[starting_column+8],splitLine[starting_column+5], splitLine[starting_column+6])
				self.addBGC(splitLine[0], splitLine[1], distance)
		
		if not self.used_ap:
			self.affinity_propagation()
		else:
			self.clusterJson()
			
	def process_network_line(self, Jaccard, DSS, AI, S, S_anchor, DDS_non_anchor, DDS_anchor):
		"""
		Processes a single line of a network file
		"""
		# Turn all strings to floats
		Jaccard = float(Jaccard)
		DSS = float(DSS)
		AI = float(AI)
		S = float(S)
		S_anchor = float(S_anchor)
		DDS_non_anchor = float(DDS_non_anchor)
		DDS_anchor = float(DDS_anchor)
		# DSS, calculated again with the new anchorboost / anchorweight
		if S_anchor != 0 and S != 0:
			non_anchor_prct = S / (S + S_anchor)
			anchor_prct = S_anchor / (S + S_anchor)

			non_anchor_weight = non_anchor_prct / (anchor_prct*self.anchor_boost + non_anchor_prct)
			anchor_weight = anchor_prct*self.anchor_boost / (anchor_prct*self.anchor_boost + non_anchor_prct)
			DDS_dis = (non_anchor_weight*DDS_non_anchor) + (anchor_weight*DDS_anchor)
		elif S_anchor == 0:
			DDS_dis = DDS_non_anchor
		else: # Only anchor domains were found
			DDS_dis = DDS_anchor

		DSS = 1.0 - DDS_dis

		# Calculate weighted values
		Jaccard_w = Jaccard*self.Jw
		DSS_w = DSS*self.DSSw
		AI_w = AI *self.AIw
		# Calculate distance between the two clusters
		distance = 1 - Jaccard_w - DSS_w - AI_w
		return distance

	def addBGC(self, gc1, gc2, distance):
		'''
		Adds a bgc to the dictionary
		'''
		self.bgcs.add(gc1)
		self.bgcs.add(gc2)
		if distance < self.cutoff:
			similarity = 1 - distance
		else:
			similarity = 0
		gcSimilarities = self.similarity_dict.setdefault(gc1,{})
		gcSimilarities[gc2] = similarity

	def affinity_propagation(self, damping=0.8):
		'''
		Does the Affinity Propagation
		'''
		# Preserve order
		self.bgcs = sorted(list(self.bgcs))
		bgc2simIdx = dict(zip(self.bgcs, range(len(self.bgcs))))
		#~ pfam_domain_categories = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Pfam-A.clans.tsv")
		#~ pfam_descrs = generatePfamDescriptionsMatrix(pfam_domain_categories)
		
		simMatrix = lil_matrix((len(bgc2simIdx), len(bgc2simIdx)), dtype=np.float32)
		for bgc1 in self.bgcs:
			# First make sure it is similar to itself
			simMatrix[bgc2simIdx[bgc1], bgc2simIdx[bgc1]] = 1
			for bgc2 in self.similarity_dict.get(bgc1,{}).keys():
				# You might get 0 values if there were matrix entries under the cutoff don't need to input these in
				# The sparse matrix
				if self.similarity_dict[bgc1][bgc2] > 1-self.cutoff:
					# Ensure symmetry
					simMatrix[bgc2simIdx[bgc1], bgc2simIdx[bgc2]] = self.similarity_dict[bgc1][bgc2]
					simMatrix[bgc2simIdx[bgc2], bgc2simIdx[bgc1]] = self.similarity_dict[bgc1][bgc2]
		
		# Use preference='min' for the minimal amount of clusters or preference='median' if you want results with more different clusters
		# In the Scikit version of AP, below, this is standard median
		labels = pysapc.SAP(damping=damping, max_iter=500, preference='min').fit_predict(simMatrix) # current
		#~ labels = pysapc.SAP(damping=0.5, max_iter=500, preference='median').fit_predict(simMatrix)
		
		# Contains the label with the corresponding clusters
		labelDict = {}
		for i in range(0, len(self.bgcs)-1):
			if labels[i] not in labelDict:
				labelDict[labels[i]] = [self.bgcs[i]]
			else:
				labelDict[labels[i]].append(self.bgcs[i])
				
		nf_name = self.ibgc.replace("\\", "/").split("/")[-1].replace(".network", "").replace(".bz2", "")
		outputFile = open("{}/affinity_propagation_{}_Jw{}_AIw{}_DSSw{}_AB{}.csv".format(self.output_dir, self.nf_name, self.Jw, self.AIw, self.DSSw, self.anchor_boost), 'w')
		for label in labelDict.keys():
			outputFile.write(str(label))
			for cluster in labelDict[label]:
				outputFile.write(","+cluster)
			outputFile.write("\n")
			
		########################################################################################################################	
		#### For json output																								####
		########################################################################################################################
		bgc_info = {}
		get_gbk_files(self.gbk_file_dir, 0, bgc_info)
		numBGCs = len(self.bgcs)
		bs_distances = [[float('%.3f' % simMatrix[row, col]) for col in xrange(row+1)] for row in xrange(numBGCs)]
		bs_data = []
		bgcJsonDict = {}
		for bgc in self.bgcs:
			bgcName = bgc
			bgcJsonDict[bgcName] = {}
			bgcJsonDict[bgcName]["id"] = bgc
			bgcJsonDict[bgcName]["desc"] = bgc_info[bgcName].description
			bgcJsonDict[bgcName]["start"] = 1
			bgcJsonDict[bgcName]["end"] = bgc_info[bgcName].max_width
			pfdFile = os.path.join(pfd_folder, bgcName + ".pfd")
			fastaFile = os.path.join(bgc_fasta_folder, bgcName + ".fasta")
			orfDict = defaultdict(dict)
			## read fasta file first to get orfs
			for line in open(fastaFile):
				if line[0] == ">":
					header = line.strip()[1:].split(':')
					if header[2]:
						orfDict[header[0]]["id"] = header[2]
					elif header[4]:
						orfDict[header[0]]["id"] = header[4]
					else:
						orfDict[header[0]]["id"] = header[0]
					orfDict[header[0]]["start"] = int(header[6])
					orfDict[header[0]]["end"] = int(header[7])
					if header[-1] == '+':
						orfDict[header[0]]["strand"] = 1
					else:
						orfDict[header[0]]["strand"] = -1
					orfDict[header[0]]["domains"] = []
			## now read pfd file to add the domains to each of the orfs
			for line in open(pfdFile):
				entry = line.split('\t')
				orf = entry[-1].strip().split(':')[0]
				pfamID = entry[5].split('.')[0]

				orfDict[orf]["domains"].append({'code': entry[5], 'start': int(entry[3]), 'end': int(entry[4]), 'bitscore': float(entry[1])})
			bgcJsonDict[bgcName]['orfs'] = orfDict.values()
		bs_data = [bgcJsonDict[bgc] for bgc in self.bgcs]
		familiesDict = {}
		for idx, label in enumerate(labels):
			members = familiesDict.setdefault(label, [])
			members.append(idx)
			familiesDict[label] = members
		bs_families = [{'id': 'FAM_%.3d' % family, 'members': members} for family, members in enumerate(familiesDict.values())]

		print("Writing JS file")
		outputFile = "{}/network_data_{}_J{}_A{}_D{}_B{}.js".format(self.output_dir, self.nf_name, self.Jw, self.AIw, self.DSSw, self.anchor_boost)
		with open(outputFile, 'w') as outfile:
			outfile.write('var bs_similarity=%s\n' % str(bs_distances))
			outfile.write('var bs_data=%s\n' % str(bs_data))
			outfile.write('var bs_families=%s' % str(bs_families))

			
######################################################################################################################
#####     #####                           SciKit affinity propagation                                  #####     #####
######################################################################################################################	

	def clusterJson(self, cutoffs=[1.0],damping=0.8):
		"""
		Another version of the affinity propagation.
		"""
		# From the data structure compute a similarity matrix for clustering, cluster using AP and then output a json
		# file with the results of the clustering for visualization
		# any distance higher than the distance cutoff will result in a similarity score of 0
		# preserve order
		self.bgcs = sorted(list(self.bgcs))
		#~ pfam_domain_categories = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Pfam-A.clans.tsv")
		#~ pfam_descrs = generatePfamDescriptionsMatrix(pfam_domain_categories)
		
		
		# because the triUdistMatrix does not accept string but only integers
		index_to_bgc = dict()
		bgc_to_index = dict()
		for i in range(len(self.bgcs)):
			bgc_to_index.update({self.bgcs[i]:i})
			index_to_bgc.update({i:self.bgcs[i]})
			
		for bgc in self.bgcs:
			if bgc in self.similarity_dict.keys():
				self.similarity_dict[bgc][bgc] = 1.0
			else:
				self.similarity_dict[bgc] = {bgc:1}
		
		for cutoff in cutoffs:
			triUdistMatrix = np.zeros((len(self.bgcs),len(self.bgcs)))
			
			for bgc1 in self.bgcs:
				for bgc2 in self.similarity_dict.get(bgc1,{}).keys():
					if self.similarity_dict[bgc1][bgc2] > 1 - cutoff:
						triUdistMatrix[bgc_to_index[bgc1]][bgc_to_index[bgc2]] = self.similarity_dict[bgc1][bgc2]
						
			symDistMatrix = triUdistMatrix + triUdistMatrix.T - np.diag(triUdistMatrix.diagonal())
			#default damping used here (=0.5)
			labels = AffinityPropagation(max_iter=500, preference=None,affinity='precomputed').fit_predict(symDistMatrix)
			familiesDict = {}
			for idx,label in enumerate(labels):
				members = familiesDict.setdefault(label,[])
				members.append(idx)
				familiesDict[label] = members

			with open("{}/new_affinity_propagation_{}_Jw{}_AIw{}_DSSw{}_AB{}.csv".format(self.output_dir, self.nf_name, self.Jw, self.AIw, self.DSSw, self.anchor_boost), 'w') as clustering_file:
				i = 0
				for label in familiesDict:
					i += 1
					clustering_file.write(str(label))
					for x in familiesDict[label]:
						clustering_file.write("," + index_to_bgc[x])
					clustering_file.write("\n")
		
		########################################################################################################################	
		#### For json output																								####
		########################################################################################################################
		#~ bgc_info = {}
		#~ get_gbk_files(self.gbk_file_dir, 0, bgc_info)
		#~ numBGCs = len(self.bgcs)
		#~ bs_distances = [[float('%.3f' % symDistMatrix[row, col]) for col in xrange(row + 1)] for row in xrange(numBGCs)]
		#~ bs_data = []
		#~ bgcJsonDict = {}
		#~ for bgc in self.bgcs:
			#~ bgcName = bgc
			#~ bgcJsonDict[bgcName] = {}
			#~ bgcJsonDict[bgcName]["id"] = bgc
			#~ bgcJsonDict[bgcName]["desc"] = bgc_info[bgcName].description
			#~ bgcJsonDict[bgcName]["start"] = 1
			#~ bgcJsonDict[bgcName]["end"] = bgc_info[bgcName].max_width
			#~ pfdFile = os.path.join(pfd_folder, bgcName + ".pfd")
			#~ fastaFile = os.path.join(bgc_fasta_folder, bgcName + ".fasta")
			#~ orfDict = defaultdict(dict)
			#~ ## read fasta file first to get orfs
			#~ for line in open(fastaFile):
				#~ if line[0] == ">":
					#~ header = line.strip()[1:].split(':')
					#~ if header[2]:
						#~ orfDict[header[0]]["id"] = header[2]
					#~ elif header[4]:
						#~ orfDict[header[0]]["id"] = header[4]
					#~ else:
						#~ orfDict[header[0]]["id"] = header[0]
					#~ orfDict[header[0]]["start"] = int(header[6])
					#~ orfDict[header[0]]["end"] = int(header[7])
					#~ if header[-1] == '+':
						#~ orfDict[header[0]]["strand"] = 1
					#~ else:
						#~ orfDict[header[0]]["strand"] = -1
					#~ orfDict[header[0]]["domains"] = []
			#~ ## now read pfd file to add the domains to each of the orfs
			#~ for line in open(pfdFile):
				#~ entry = line.split('\t')
				#~ orf = entry[-1].strip().split(':')[0]
				#~ pfamID = entry[5].split('.')[0]
				#~ pfamDescr = pfam_descrs.get(pfamID,None)
				#~ if pfamDescr:
					#~ orfDict[orf]["domains"].append({'code': '{} : {}'.format(pfamID,pfamDescr),'start':int(entry[3]),'end':int(entry[4]),'bitscore': float(entry[1])})
				#~ else:
				#~ orfDict[orf]["domains"].append({'code': entry[5], 'start': int(entry[3]), 'end': int(entry[4]), 'bitscore': float(entry[1])})
			#~ bgcJsonDict[bgcName]['orfs'] = orfDict.values()
		#~ bs_data = [bgcJsonDict[bgc] for bgc in self.bgcs]
		#~ familiesDict = {}
		#~ for idx, label in enumerate(labels):
			#~ members = familiesDict.setdefault(label, [])
			#~ members.append(idx)
			#~ familiesDict[label] = members
		#~ bs_families = [{'id': 'FAM_%.3d' % family, 'members': members} for family, members in enumerate(familiesDict.values())]

		#~ print("  Writing JS file")
		#~ outputFile = "{}/network_data_{}_J{}_A{}_D{}_B{}.js".format(self.output_dir, self.nf_name, self.Jw, self.AIw, self.DSSw, self.anchor_boost)
		#~ with open(outputFile, 'w') as outfile:
			#~ outfile.write('var bs_similarity=%s\n' % str(bs_distances))
			#~ outfile.write('var bs_data=%s\n' % str(bs_data))
			#~ outfile.write('var bs_families=%s' % str(bs_families))
		


if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='Affinity_Popagation_new_weights.py')
	# To be sure nothing will be cut off when no cutoff is given, a default is set to 1.0
	parser.add_argument("-c" ,"--cut", 
	help="the cutoff values seperated by a komma (default = 1.0)", required=False, default="1.0")
	parser.add_argument("-n", "--network", 
	help="The BGC network file", required=True)
	parser.add_argument("-j", "--jaccard_weight", 
	help="The weight for the jaccard index", required=True)
	parser.add_argument("-a", "--adjacency_index_weight", 
	help="The weight for the adjacency index", required=True)
	parser.add_argument("-d", "--domain_similarity_score_weight", 
	help="The weight for the domain similarity score", required=True)
	parser.add_argument("-b", "--anchor_boost", 
	help="The weight for the anchor boost", required=True)
	parser.add_argument("-o", "--output_dir", 
	help="Directory to put the output file(s) in. (default ./)", required=False, default="./")
	parser.add_argument("-p", "--new_affinity_propagation", 
	help="If this flag is used the scikit affinity propagation will be used",action="store_true", default=False)
	parser.add_argument("-g", "--gbk_files", 
	help="The directory with gbk_files", default="/home/roete009/Thesis/python/data_files/genbank_files", required=False)
	parser.add_argument("-f", "--data_folder", 
	help="The directory fasta, pfd, psd, etc.", default="/mnt/scratch/navar016/Actinobacteria_analyzed/BiG-SCAPE_results/2017-03-02_Actinobacteria")
	
	args = parser.parse_args()
	try:
		# Checks for the cutoff values
		if "," in args.cut:
			cutoffs = args.cut.split(",")
		else:
			cutoffs = [args.cut]
		# And if all the given cutoffs are valid numbers
		for cutoff in cutoffs:
			float(cutoff)
	except:
		sys.exit("The given parameters cutoff values are not correct!")
		
	try:
		if float(args.jaccard_weight) + float(args.adjacency_index_weight) + float(args.domain_similarity_score_weight) != float(1.0):
			print("Warning: your weights did not add up to 1.0")
		
		if float(args.anchor_boost) < 1.0:
			print("Warning: the anchor boost is smaller than 1")
	except TypeError as e:
		sys.exit("You did not provide a valid number for one of the weights: " + e)
		
	
	try:
		if os.path.isdir(args.data_folder):
			bgc_fasta_folder = os.path.join(args.data_folder, "fasta")
			pfd_folder = os.path.join(args.data_folder, "pfd")
		else:
			sys.exit("provide a folder with all fasta files")
		if not os.path.isdir(args.output_dir):
			sys.exit("specified output directory " + args.output_dir + " is not valid")
		# Check if it's a network file and if it can be read
		if ".network" in args.network or ".bz2" in args.network:
			if os.path.isdir(args.gbk_files):
				for cutoff in cutoffs:
					# InputFile, jaccard_weight, adjacency index weight, domain similarity score weight, anchorboost, cutoff
					NFR = NetworkFileReader(args.network, float(args.jaccard_weight), float(args.adjacency_index_weight), 
					float(args.domain_similarity_score_weight), float(args.anchor_boost), args.output_dir, cutoff, args.new_affinity_propagation, args.gbk_files)
					NFR.process_network_file()
			else:
				sys.exit("Please provide a folder with gbk files of the BGCs")
		else:
			sys.exit("The file you provided has to be a network file")
			# Check if valid output directory
	except IOError:
		sys.exit("The network file or the output directory is not valid!\n" + args.network)












