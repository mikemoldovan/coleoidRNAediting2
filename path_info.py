"""
path_info
"""

import re
from os import listdir, system, mkdir, path
from Bio import AlignIO, SeqIO
from ete3 import Tree
from optparse import OptionParser

#treetime ancestral --aln oct_als/orfs_oct_comp130875_c0_seq10.fasta --tree trees/orfs_oct_comp130875_c0_seq10.fasta.phy.treefile --outdir tmp/

#We have output directory, working with it
def metadata_parse(meta_file):
    meta_dict = {}
    with open(meta_file) as inh:
        for s in inh:
            s = s.strip().split()
            s_name = s[0] + s[2]
            s[1] = list(map(eval, s[1].split(',')))
            meta_dict[s_name] = s
    return meta_dict


def get_seqs(fafile):
    fa_dict = {}
    with open(fafile) as inh:
        for rec in SeqIO.parse(fafile, "fasta"):
            s = ""
            for c in rec.seq:
                s += c
            fa_dict[rec.id] = s
    return fa_dict

def get_tree(treefile):
    with open(treefile) as inh:
        skip = True
        for s in inh:
            if s.startswith("Begin Trees;"):
                skip = False
            if not skip and s.startswith("End;"):
                skip=True
            if skip:
                continue
            
            if not "tree1=" in s:
                continue

            s = s.split("tree1=")
            treestr = ""
            is_G2A = True if """&mutations="G""" in s[1] else False
            
            right_format_char = True
            for c in s[1]:
                if c == "[":
                    right_format_char = False
                if right_format_char:
                    treestr += c
                if c == "]":
                    right_format_char = True
    t = Tree(treestr, format=1)
    return t, is_G2A


def differences(node, par_node, seqs, clust_data):
    s1 = seqs[node.name]
    s2 = seqs[par_node.name]
    diff_list = []
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            diff_list.append(clust_data[1][i])
    return diff_list
    

def process_tree(seqs, tree, clust_data, outh):
    for leaf in tree:
        node = leaf
        subst_list = []
        while not node.is_root():
            par_node = node.up
            diff_list = ','.join(map(str, differences(node, par_node, seqs, clust_data)))
            subst_list.append(diff_list)
            node = node.up
        
        outh.write("{}\t{}\t{}\t{}\n".format(clust_data[0], clust_data[2], leaf.name.split('_')[-1], '|'.join(subst_list)))
            


def mainloop(meta_file, als_dir, trees_dir, out_dir, outfile, treetime_run = True):
    meta_dict = metadata_parse(meta_file)
    
    with open(outfile, 'w') as outhandle:
        for clust_id, clust_data in meta_dict.items():
            al_filename = als_dir + clust_id + ".fasta"
            tree_filename = trees_dir + clust_id + ".fasta.phy.treefile"
            if not path.isfile(tree_filename):
                continue
    
            outdir_name = out_dir + clust_id + '/'
            
            if treetime_run:
                system("treetime ancestral --aln {} --tree {} --outdir {}".format(al_filename, tree_filename, outdir_name))
            
            seqs = get_seqs(outdir_name + "ancestral_sequences.fasta")
            tree, is_G2A = get_tree(outdir_name + "annotated_tree.nexus")
            
            if is_G2A:
                continue
            
            process_tree(seqs, tree, clust_data, outhandle)


parser = OptionParser()
parser.add_option("-m", "--meta_file", help="File with cluster metadata")
parser.add_option("-a", "--als_dir", help="Directory with alignments")
parser.add_option("-t", "--trees_dir", help="Directory with trees")
parser.add_option("-n", "--out_dir", help="Directory for/with treetime output")
parser.add_option("-o", "--outfile", help="Output file")
parser.add_option("-r", "--treetime_run", help="Need a treetime run (Y or N)?. def. Y", default='Y')
opt, args = parser.parse_args()

mainloop(opt.meta_file, 
         opt.als_dir, 
         opt.trees_dir, 
         opt.out_dir, 
         opt.outfile, 
         True if opt.treetime_run == 'Y' else False) 




