"""
call_edsite_variants

For each read mapped onto the editing site, determine the corresponding nucleotide in the read.
"""

from optparse import OptionParser


def readfasta(fafile):
    fastadict = dict()
    seq = ""
    with open(fafile) as inhandle:
        for s in inhandle:
            if s[0] == '>':
                if seq:
                    fastadict[seq_id] = seq
                    seq = ""
                seq_id = s.strip().split()[0][1:]
            else:
                seq += s.strip()
    fastadict[seq_id] = seq
    return fastadict


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


def parse_cigar(cigar_str):
    cigar_arr = []
    elem = ""
    for c in cigar_str:
        if c.isdigit():
            elem += c
        else:
            elem += c
            cigar_arr.append(elem)
            elem = ""
    return cigar_arr

def format_read(read_seq, cigar_arr):
    crd = 0
    for c in cigar_arr:
        c_com = c[-1]
        c_len = eval(c[:-1])
        if c_com == "M":
            crd += c_len
        elif c_com == "D":
            read_seq = read_seq[:crd] + "N"*c_len + read_seq[crd:]
            crd += c_len
        elif c_com == "I":
            read_seq = read_seq[:crd] + read_seq[crd + c_len:]
        elif c_com == "S":
            read_seq = read_seq[:crd] + read_seq[crd + c_len:]
        else:
            print(c_com)
    return read_seq

def edsite_cov(ref_seq_id, 
               ref_seq, 
               read_id, 
               read_seq, 
               cigar_arr, 
               align_start, 
               edsite_crds):
    
    good_al = False
    for elem in cigar_arr:
        if elem.endswith("M"):
            if eval(elem[:-1]) > 80:
                good_al = True
    if not good_al:
        return 0
    
    new_read_seq = format_read(read_seq, cigar_arr)
    l = len(new_read_seq)
    for crd in edsite_crds.keys():
        crd_on_read = crd - align_start
        if crd_on_read >= l:
            continue
        print("{}\t{}\t{}\t{}\t{}".format(ref_seq_id, read_id, crd, ref_seq[crd - 1], new_read_seq[crd_on_read]))
#    print(cigar_arr, ref_seq_id, read_id)
#    ref_alseq = ref_seq[align_start-1:align_start + l-1]
#    print(new_read_seq)
#    print(ref_alseq)
#    print("")
    
#    align_crds = []
#    print(read_id, cigar_arr)
#    return 0

def edsites_cov(sam_file, transcr_fasta, edsite_crds):
    edsite_crd_dict = read_edsite_crd_file(edsite_crds, syn_nsyn = None) 
    fastadict = readfasta(transcr_fasta)
    with open(sam_file) as inhandle:
        for s in inhandle:
            s = s.strip().split()
            ref_seq_id = s[2]
            if not edsite_crd_dict.get(ref_seq_id):
                continue
            ref_seq = fastadict[ref_seq_id]
            read_id = s[0]
            read_seq = s[9]
            cigar_arr = parse_cigar(s[5])
            align_start = eval(s[3])
            seq_edsites = edsite_crd_dict[ref_seq_id]
            
            edsite_crds = dict()
            for edsite_crd in seq_edsites.keys():
                if edsite_crd >= align_start - 1 and edsite_crd <= align_start + 160:
                    edsite_crds[edsite_crd] = True
            if not edsite_crds:
                continue
                
            edsite_cov(ref_seq_id, 
                       ref_seq, 
                       read_id, 
                       read_seq, 
                       cigar_arr, 
                       align_start, 
                       edsite_crds)


parser = OptionParser()
parser.add_option("-s", "--sam_file", help="SAM file with the alignment")
parser.add_option("-t", "--transcr_fasta", help="Corresponding FASTA-file with transcripts")
parser.add_option("-e", "--edsite_crds", help="Editing site coordinates")
opt, args = parser.parse_args()

edsites_cov(opt.sam_file, opt.transcr_fasta, opt.edsite_crds)