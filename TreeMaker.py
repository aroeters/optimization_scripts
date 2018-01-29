import argparse
import sys
import os
import random

# Imports for creating the phylogenetic Tree
from PyQt4 import QtCore
from PyQt4.QtGui import QGraphicsRectItem, QColor, QBrush, QGraphicsTextItem, QGraphicsEllipseItem, QPen, QLabel
from ete3 import Tree, TreeStyle, faces
import ete3

#-i /mnt/scratch/roete009/big_scape_output/networks_all_hybrids_lcs -c 0.50 -o ./ -p ./genome_plasmids_matches.txt -n ./complete_Streptomyces_tree.nwk

# A lambda function to get a random number between 0 and 255 to be used for the rgb color for the nodes
r = lambda: random.randint(0,255)
# Extend the given TreeNode class to be able to use more attributes
class ExtendedTreeNode(ete3.coretype.tree.TreeNode):
    def add_info(self, new_info):
        try:
            self.info += "," + new_info
        except:
            self.info = new_info

    def get_info(self):
        try:
            return self.info
        except:
            return ""

    def set_class(self, bgc_class):
        self.bgc_class = bgc_class

    def get_class(self):
        try:
            return self.bgc_class
        except:
            return ""

    def set_visible_classes(self, visible):
        self.visible_classes = visible

    def get_visible_classes(self):
        try:
            return self.visible_classes
        except:
            return set()

    def set_color_code(self, color_code):
        self.color_code = color_code

    def get_color_code(self):
        try:
            return self.color_code
        except:
            return "#FFFFFF" # white

class TreeMaker():
    """
    This class will make a phylogenetic tree with the biosynthetic gene clusters (BGCs) next to the tree.
    By doing so, missing or extra BGCs per species can be detected and horizontal BGC transfer might be seen.
    """
    def __init__(self):
        self.output_name = "heatmap_input_c80_NRPS.csv"
        self.classes_to_show = {"NRPS","PKSI","PKS-NRP"}
        self.get_cmd_arguments()
        self.create_phylo_tree()
        self.get_bgc_cluster_families()
        self.match_node_names_to_accession()
        self.get_plasmids_accessions_per_genome_accession()
        self.add_info_to_node()
        self.add_faces_to_nodes()
        
    def get_cmd_arguments(self):
        """
        Retrieves the command line input for this script to run
        """
        print("Retrieving command line arguments...")
        parser = argparse.ArgumentParser(prog='network_connectivity.py')
        parser.add_argument("-i", "--big_scape_output",
        help="The directory that contains all network files of BiG-SCAPE\nWill iterate over sub-folders to find the needed files.\nThis should end with /networks_all or networks_all_hybrids_lcs", required=True)
        
        parser.add_argument("-c", "--cutoff",
        help="The cutoff value that should be used on two BGCs as a maximum distance\nfor which one would still call them the same BGC.", required = False, default="0.5")
        
        parser.add_argument("-o", "--output_dir",
        help="The directory to output the final results in.", required=False, default="./")
        
        parser.add_argument("-p", "--plasmids",
        help="The file that contains the organisms accession with the plasmids accessions\n(if there are any) in a tab delimited file.", required=True)
        
        parser.add_argument("-n", "--newick",
        help="The file that contains phylogenetic tree.", required=True)
        
        parser.add_argument("-m", "--with_mibig",
        help="Only show the cluster families that also have an MIBiG cluster.", required=False, action='store_true')
        args = parser.parse_args()

        if not os.path.isfile(args.plasmids):
            sys.exit("You have not provided a valid file for the plasmids:\n{}".format(args.plasmids))
        if not os.path.isfile(args.newick):
            sys.exit("Please provide a newick formatted file containing the wanted phylogenetic tree:\n{}".format(args.newick))
        if not os.path.isdir(args.output_dir):
            sys.exit("Please provide a valid output directory:\n{}".format(args.output_dir))
        if not os.path.isdir(args.big_scape_output):
            sys.exit("Please provide a valid directory for the BiG-SCAPE output:\n{}".format(args.big_scape_output))
        try:
            if args.cutoff.endswith("0") or args.cutoff.endswith("5"):
                if float(args.cutoff) >= 0.20 and float(args.cutoff) <= 0.80:
                    self.big_scape_files = self.get_cluster_files(args.big_scape_output, args.cutoff)
                else:
                    sys.exit("Provide a number between 0.20 and 0.80")
            else:
                sys.exit("Provide a number between 0.20 and 0.80 that is a multiple of 0.05")
        except ValueError:
            sys.exit("Please provide a valid number bewteen 0.20 and 0.80 that is a multiple of 0.05")

        self.outdir = args.output_dir
        self.plasmids_file = args.plasmids
        self.nwk = args.newick
        self.mibig_related = args.with_mibig

    def get_cluster_files(self, big_scape_dir, cutoff):
        """
        Finds all the needed cluster files that shows which BGCs belong together
        """
        list_of_files = {}
        for (dirpath, dirnames, filenames) in os.walk(big_scape_dir):
            for filename in filenames:
                if filename.endswith('.tsv') and cutoff in filename:
                    bgc_class = filename.split("_")[0].replace("all","")
                    list_of_files[bgc_class] = os.sep.join([dirpath, filename])
        return list_of_files

    def create_phylo_tree(self):
        """
        Creates a basic tree with only the leafs and branches by reading the tree file in nwk format.
        """
        self.ts = TreeStyle()
        self.ts.show_leaf_name = True
        self.ts.layout_fn = master_ly
        self.ts.title.add_face(faces.TextFace("Phyolgenetic tree showing BGCs per organism", fsize=25), 0)
        self.tree = Tree(self.nwk)
    
    def get_bgc_cluster_families(self):
        """
        Parses all the needed files coming from BiG-SCAPE and retrieves all the clusters it has made.
        These clusters of bgcs should in theory be the same (or atleast extremely similar) bgcs in another organism.
        """

        self.cluster_families = dict()
        for bgc_class in self.big_scape_files:
            if bgc_class in self.classes_to_show:
                file_in = self.big_scape_files[bgc_class]
                self.cluster_families.update({bgc_class:dict()})

                with open(file_in) as opened_file:
                    opened_file.readline() # to skip the header in the files
                    for line in opened_file:
                        BGC, family = line.strip().split()
                        if family in self.cluster_families[bgc_class]:
                            self.cluster_families[bgc_class][family].append(BGC)
                        else:
                            self.cluster_families[bgc_class].update({family:[BGC]})
                            
    def match_node_names_to_accession(self):
        """
        Matches the accession name to the node name
        """
        self.accession_to_node_name = dict()
        for node in self.tree:
            if node.is_leaf():
                accession = node.name.split("_")[-1]
                self.accession_to_node_name.update({accession:node.name})

    def get_plasmids_accessions_per_genome_accession(self):
        """
        Reads in the plasmids file and extracts the plasmids per genome.
        """
        for line in open(self.plasmids_file):
            split_line = line.strip().split("\t")
            genome_accession = split_line[0]
            if len(split_line) >= 2:
                if genome_accession in self.accession_to_node_name:
                    node_name = self.accession_to_node_name[genome_accession]
                else:
                    node_name = "-"
                for plasmid_accession in split_line[1::]:
                    self.accession_to_node_name.update({plasmid_accession:node_name})

    def add_info_to_node(self):
        """
        Adds the info of all the clusters to the names of the species for later handling
        """
        used_colors = set()
        absence_presence = {self.accession_to_node_name[accession] : dict() for accession in self.accession_to_node_name.keys()}
        bgcs = set()
        for bgc_class in self.cluster_families:
            for cluster_family in self.cluster_families[bgc_class]:
                cluster_bgcs = self.cluster_families[bgc_class][cluster_family]

                # Filter families that only contain MIBiG BGCs and no clusters from the organisms
                # The families left only contain a BGC from MIBiG with other from an organism or 
                # only from the organisms that have been put through antiSMASH 
                mibig_bgcs = [x.split(".")[0] for x in cluster_bgcs if "BGC" in x]
                if self.mibig_related and len(mibig_bgcs) != 0 and len(mibig_bgcs) != len(cluster_bgcs):
                    cluster_family_color = '#%02X%02X%02X' % (r(),r(),r())
                    while cluster_family_color in used_colors:
                        cluster_family_color = '#%02X%02X%02X' % (r(),r(),r())
                    used_colors.add(cluster_family_color)
                    mibig_names = "-"
                    if len(mibig_bgcs) != 0:
                        mibig_names = ";".join(mibig_bgcs)
                    for bgc in cluster_bgcs:
                        self.update_nodes(bgc, bgc_class, mibig_names, cluster_family_color, cluster_family)
                    
                    if self.mibig_related:
                        absence_presence_name = "_".join(mibig_bgcs)
                        for accession in self.accession_to_node_name.keys():
                            accessions_in_cluster = [x.split(".")[0] for x in cluster_bgcs]
                            if accession in accessions_in_cluster:
                                absence_presence[self.accession_to_node_name[accession]].update({absence_presence_name:1})
                            else:
                                try:
                                    absence_presence[self.accession_to_node_name[accession]][absence_presence_name]
                                except:
                                    absence_presence[self.accession_to_node_name[accession]].update({absence_presence_name:0})
                            bgcs.add(absence_presence_name)

        self.write_heatmap_input(absence_presence, sorted(list(bgcs)))

    def write_heatmap_input(self, absence_presence, bgcs):
        """
        Writes the input for a heatmap version of the phylogenetic tree
        """
        heatmap_file = open(self.output_name, 'w')
        heatmap_file.write("species," + ",".join(bgcs) + "\n")
        for species in absence_presence:
            heatmap_file.write(species)
            bgc_string = ""
            for bgc in bgcs:
                bgc_string += ",{}".format(absence_presence[species][bgc])
            heatmap_file.write(bgc_string + "\n")
        
        
    def update_nodes(self, bgc, bgc_class, mibig_bgcs, cluster_family_color, cluster_family):
        """
        Updates the nodes with the info they need, used by add_info_to_node()
        """
        if not bgc.startswith("BGC"):
            accession, version, cluster = bgc.split(".")
            if accession in self.accession_to_node_name:
                node_info = self.accession_to_node_name[accession]
                node_list = self.tree.search_nodes(name=node_info)
                if len(node_list) != 0: # A node is found
                    node = node_list[0]
                    node.__class__ = ExtendedTreeNode
                    node.add_info("{}:{}:{}:{}:{}".format(bgc_class, cluster_family, mibig_bgcs, cluster, cluster_family_color))
                    node.set_class(bgc_class)
                    node.set_visible_classes(self.classes_to_show)
    
    def add_faces_to_nodes(self):
        """
        Adds the faces to the corresponding nodes
        """
        self.ts.layout_fn = master_ly

    def show_tree(self):
        """
        Shows the tree with the needed style
        """
        self.tree.show(tree_style=self.ts)

class InteractiveItem(QGraphicsRectItem):
    def __init__(self, bgc_info, *args, **kargs):
        QGraphicsRectItem.__init__(self, *args, **kargs)
        self.node = None
        self.label = None
        self.bgc_info = bgc_info
        self.setCursor(QtCore.Qt.PointingHandCursor)
        self.setAcceptsHoverEvents(True)

    def hoverEnterEvent (self, e):
        # Show/hide a text item over the custom DynamicItemFace
        if not self.label:
            self.label = QLabel(self.bgc_info)
        self.label.show()

    def hoverLeaveEvent(self, e):
        if self.label:
            self.label.hide()

def bgc_name_face(node, *args, **kargs):
    """
    This is the item generator. It must receive a node object, and
    returns a Qt4 graphics item that can be used as a node face.
    """
    # Receive an arbitrary number of arguments, in this case width and
    # Height of the faces and the information about the BGC
    width = args[0]
    height = args[1]
    # Add the popup
    interactive_face = InteractiveItem("Class : {}\nRelated MIBiG : {}\nCluster family : {}".format(args[2], args[4], args[3]), 0, 0, width, height)
    # Keep a link within the item to access node info
    interactive_face.node = node
    # Remove border around the masterItem
    interactive_face.setPen(QPen(QtCore.Qt.NoPen))
    # Add ellipse around text
    ellipse = QGraphicsEllipseItem(interactive_face.rect())
    ellipse.setParentItem(interactive_face)
    # Change ellipse color
    ellipse.setBrush(QBrush(QColor(args[6])))
    # Add node name within the ellipse
    text = QGraphicsTextItem(args[5])
    text.setTextWidth(50)
    text.setParentItem(ellipse)
    # Center text according to masterItem size
    text_width = text.boundingRect().width()
    text_height = text.boundingRect().height()
    center = interactive_face.boundingRect().center()
    text.setPos(center.x()-text_width/2, center.y()-text_height/2)
    return interactive_face

def master_ly(node):
    # Change the class of the node to the custom made class with extra functionality
    node.__class__ = ExtendedTreeNode
    if node.is_leaf():
        # Adds the face to the node in the Tree
        try:
            for i in range(0,node.info.count(",") + 1):
                bgc_class, cluster_family, mibig_clusters, cluster, color = node.info.split(",")[i].split(":")
                face = faces.DynamicItemFace(bgc_name_face, 100, 35, bgc_class, cluster_family, mibig_clusters, cluster, color)
                node.add_face(face, position="aligned", column = i)
        except:
            # if the node has no bgcs that are linked to a mibig family (-m mode) or there are no bgcs at all
            # add no faces to that node
            pass

def main():
    treemaker = TreeMaker()
    treemaker.show_tree()

if __name__ == "__main__":
    main()
