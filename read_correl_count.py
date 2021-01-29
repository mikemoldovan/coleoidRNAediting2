"""
read_correl_count

Calculate correlations in editing site data based on the read mapping
Works with call_edsite_variants.py  output

Input:
1. Editing sites table
2. call_edsite_variants.py output table
"""

from optparse import OptionParser

class Correl_Pos():
	def __init__(self, dist):
		self.AA = 0
		self.AG = 0
		self.GA = 0
		self.GG = 0
		self.dist = dist


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


def make_pair_dict(edsite_crd_dict, radius=200):
	pair_dict = dict()
	for seq_id in edsite_crd_dict.keys():
		pair_dict[seq_id] = dict()
		c_arr = sorted(edsite_crd_dict[seq_id].keys())
		for crd1 in c_arr:
			for crd2 in c_arr:
				if crd2 <= crd1:
					continue
				dist = crd2 - crd1
				if dist > radius:
					continue
				pair_dict[seq_id][(crd1,crd2)] = Correl_Pos(dist)
	return pair_dict


def count_values(oneread_lets, pair_dict):
	crds = sorted(oneread_lets.keys())
	l = len(crds)
	if l < 2:
		return 0
	for i in range(l):
		for j in range(i+1,l):
			pair = (crds[i],crds[j])
#			if not pair_dict.get(pair):
#				print(pair_dict, pair)
			if oneread_lets[crds[i]] == 'A' and oneread_lets[crds[j]] == 'A':
				pair_dict[pair].AA += 1
			elif oneread_lets[crds[i]] == 'A' and oneread_lets[crds[j]] == 'G':
				pair_dict[pair].AG += 1
			elif oneread_lets[crds[i]] == 'G' and oneread_lets[crds[j]] == 'A':
				pair_dict[pair].GA += 1
			elif oneread_lets[crds[i]] == 'G' and oneread_lets[crds[j]] == 'G':
				pair_dict[pair].GG += 1


def parse_mism_calls(mism_calls_tsv, pair_dict):
	with open(mism_calls_tsv) as inhandle:
		oneread_lets = dict()
		curr_read = ""
		curr_gene = ""
		count = 0
		for s in inhandle:
			count += 1
#			if count %10000 == 0:
#				print(count)
			s = s.strip().split()
			if curr_read != s[1]:
				pairs = pair_dict.get(curr_gene)
				if pairs:
					count_values(oneread_lets, pairs)
				curr_read = s[1]
				curr_gene = s[0]
				oneread_lets = dict()
			oneread_lets[eval(s[2])] = s[4]


def print_res(pair_dict):
	for seq_id in pair_dict.keys():
		for (crd1, crd2) in pair_dict[seq_id].keys():
			obj = pair_dict[seq_id][(crd1, crd2)]
			cov = float(obj.AA*obj.GG) - float(obj.AG*obj.GA)
			norm = float(obj.AA + obj.AG)*(obj.AA + obj.GA)*(obj.GG + obj.GA)*(obj.GG + obj.AG)
			if norm == 0:
				continue
			print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(seq_id, 
																	  crd1, 
																	  crd2, 
																	  obj.dist, 
																	  obj.AA, 
																	  obj.AG, 
																	  obj.GA, 
																	  obj.GG, 
																	  cov, 
																	  norm, 
																	  cov/(norm**0.5)))

def main(edsite_crd_file, mism_calls_tsv):
	edsite_crd_dict = read_edsite_crd_file(edsite_crd_file, syn_nsyn = None)
	pair_dict = make_pair_dict(edsite_crd_dict, radius=200)
	parse_mism_calls(mism_calls_tsv, pair_dict)
	print_res(pair_dict)


parser = OptionParser()
parser.add_option("-e", "--edsite_crd_file", help="Editing sites table")
parser.add_option("-m", "--mism_calls_tsv", help="call_edsite_variants.py output table")
opt, args = parser.parse_args()


main(opt.edsite_crd_file, opt.mism_calls_tsv)
