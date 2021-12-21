"""
build_state_als

Builds alignments of editing states based on read mapping data
"""


import numpy as np
from optparse import OptionParser

def read_edsite_crd_file(edsite_crd_file, syn_nsyn = None):
    edsite_crd_dict = dict()
    with open(edsite_crd_file) as inhandle:
        for s in inhandle:
            if s[0] == '#':
                continue
            if syn_nsyn == "syn" and s[4] != "syn":
                continue
            elif syn_nsyn == "nsyn" and s[4] == "syn":
                continue
            s = s.strip().split()
            if not edsite_crd_dict.get(s[0]):
                edsite_crd_dict[s[0]] = dict()
            edsite_crd_dict[s[0]][eval(s[3])] = s
    return edsite_crd_dict


def window_search(edsite_crd_dict, window_size = 50, min_edsite_num = 4):
    windows = {}
    for seq_id, crd_dict in edsite_crd_dict.items():
        crd_list = sorted(crd_dict.keys())
        l = len(crd_list)
        lookup_list = [False for i in range(l)]
        found_site = True
        while found_site:
            longest_clust_len = 0
            longest_clust_crds = None
            found_site = False
            for i in range(l - min_edsite_num):
                if crd_list[i + min_edsite_num] - crd_list[i] > window_size:
                    continue
                j = min_edsite_num
                while i + j < l and crd_list[i + j] - crd_list[i] <= window_size:
                    j += 1
                if np.sum(lookup_list[i:i + j]) > 0:
                    continue
                clust = crd_list[i:i+j]
                clust_len = j
                if clust_len > longest_clust_len:
                    longest_clust_len = clust_len
                    longest_clust_crds = i, i+j
                    found_site = True
            
            if found_site:
                if not windows.get(seq_id):
                    windows[seq_id] = []
                windows[seq_id].append(crd_list[longest_clust_crds[0]:longest_clust_crds[1]])
                for i in range(longest_clust_crds[0], longest_clust_crds[1]):
                    lookup_list[i] = True

    return windows


def read_info(align_tab):
    read_id = None
    with open(align_tab) as inh:
        for s in inh:
            s = s.strip().split()
            if not read_id:
                read_id = '.'.join(s[1].split('.')[:3])
                seq_id = s[0]
                covered_nucls = {}
                good_read = True
            if '.'.join(s[1].split('.')[:3]) != read_id:
                yield read_id, seq_id, covered_nucls, good_read
                read_id = '.'.join(s[1].split('.')[:3])
                seq_id = s[0]
                covered_nucls = {}
                good_read = True
            covered_nucls[eval(s[2])] = s[4]
            if s[4] not in "AG":
                good_read = False
    yield read_id, seq_id, covered_nucls, good_read

    
def make_read_seq(clust, covered_nucls):
    seq = ""
    for i in clust:
        if not covered_nucls.get(i):
            return None, None
        seq += covered_nucls[i]
    return seq, tuple(clust)

def filter_clusters_by_coverage(align_tab, windows_dict):
    windows_cov_dict = {} #seq_id -> window_pos_tuple -> seq -> num_seq
    for read_id, seq_id, covered_nucls, good_read in read_info(align_tab):
        if not good_read:
            continue
        wind_list = windows_dict.get(seq_id)
        if not wind_list:
            continue
        read_in_clust = False
        minval = np.min(list(covered_nucls.keys()))
        maxval = np.max(list(covered_nucls.keys()))
#        print(covered_nucls, minval, maxval, read_id)
        clust_list = []
        for clust in wind_list:
            if minval <= clust[0] and maxval >= clust[1]:
                read_in_clust = True
                clust_list.append(clust)
        if not read_in_clust:
            continue
        
        if not windows_cov_dict.get(seq_id):
            windows_cov_dict[seq_id] = dict()
        
        for clust in clust_list:
            read_seq, t_clust = make_read_seq(clust, covered_nucls)
            if not read_seq:
                continue
            if not windows_cov_dict[seq_id].get(t_clust):
                windows_cov_dict[seq_id][t_clust] = {}
            if not windows_cov_dict[seq_id][t_clust].get(read_seq):
                windows_cov_dict[seq_id][t_clust][read_seq] = [0, read_id]
            windows_cov_dict[seq_id][t_clust][read_seq][0] += 1
    return windows_cov_dict
            
#Add all-As
#Filter by support of each state -- 5 reads
#Filter by state number -- minimum 3 without outgroup

def print_alignment(alignment_dir, read_seqs, count, seq_id, readnum_cutoff):
    al_dict = {}
    was_ref = False
    for seq, info in read_seqs.items():
        l = len(seq)
        if info[0] < readnum_cutoff:
            continue
        if not 'G' in seq:
            was_ref = True
            info[1] = "out_" + info[1]
        name = info[1] + '_' + str(info[0])
        al_dict[name] = seq
        
    if not was_ref:
        al_dict["out"] = 'A'*l
    
    if len(al_dict.keys()) < 3:
        return None
    
    with open(alignment_dir + seq_id + str(count) + ".fasta", 'w') as outh:
        for name, seq in al_dict.items():
            outh.write('>' + name + '\n')
            outh.write(seq + '\n')     

            
def make_alignments(windows_cov_dict, metadata_filename, alignment_dir, readnum_cutoff=5):
    with open(metadata_filename, 'w') as outh:
        for seq_id, seq_info in windows_cov_dict.items():
            count = 0
            for clust_crds, read_seqs in seq_info.items():
                outh.write(seq_id + '\t' + ','.join(map(str, list(clust_crds))) + '\t' + str(count) + '\n')
                print_alignment(alignment_dir, read_seqs, count, seq_id, readnum_cutoff)
                count += 1

                
def prepare_data(edsite_crd_file, 
                 align_tab, 
                 metadata_filename, 
                 alignment_dir, 
                 readnum_cutoff=5, 
                 window_size=50, 
                 min_edsite_num=4):

    edsite_crd_dict =  read_edsite_crd_file(edsite_crd_file)   
    windows_dict = window_search(edsite_crd_dict, window_size, min_edsite_num)
#    windows_dict = {'orfs_oct_comp183909_c0_seq7': [[8663, 8664, 8671, 8674, 8675, 8685, 8694, 8695, 8704], [10598, 10600, 10601, 10608, 10620, 10621, 10642, 10643, 10644], [11625, 11633, 11637, 11639, 11646, 11666, 11667, 11671], [11073, 11075, 11092, 11093, 11110, 11111, 11119], [14358, 14361, 14370, 14377, 14390, 14392, 14402], [7000, 7006, 7007, 7021, 7033, 7036], [7818, 7835, 7836, 7839, 7844, 7845], [10239, 10246, 10247, 10252, 10255, 10274], [11416, 11419, 11429, 11436, 11438, 11456], [15174, 15187, 15210, 15212, 15213, 15214], [6741, 6756, 6783, 6785, 6789], [9058, 9059, 9062, 9071, 9072], [10142, 10145, 10149, 10154, 10157], [10461, 10466, 10484, 10500, 10503], [10875, 10878, 10899, 10906, 10922], [11253, 11258, 11289, 11293, 11295], [12912, 12916, 12942, 12956, 12957], [13253, 13254, 13278, 13281, 13299]]}

    windows_cov_dict = filter_clusters_by_coverage(align_tab, windows_dict)
    make_alignments(windows_cov_dict, metadata_filename, alignment_dir, readnum_cutoff=5)



parser = OptionParser()
parser.add_option("-p", "--edsite_crd_file", help="File with editing site coordinates")
parser.add_option("-l", "--align_tab", help="Table with read mapping output")
parser.add_option("-i", "--metadata_filename", help="Name of the generated file with metadata")
parser.add_option("-a", "--alignment_dir", help="Directory with generated alignments")
parser.add_option("-b", "--readnum_cutoff", help="Minimal number of state-supporting reads. def. 5", default='5')
parser.add_option("-c", "--window_size", help="Maximal length of a cluster. def. 50nt", default='50')
parser.add_option("-d", "--min_edsite_num", help="Minimal number of editing sites per cluster -1 . def. 4 (for 5 sites)", default='4')
opt, args = parser.parse_args()



prepare_data(opt.edsite_crd_file, 
             opt.align_tab, 
             opt.metadata_filename, 
             opt.alignment_dir, 
             eval(opt.readnum_cutoff), 
             eval(opt.window_size), 
             eval(opt.min_edsite_num))

