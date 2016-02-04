"""
Implementation of Feng Zhang's scoring algorithm for mismatches between guideRNAs and off-target sites,
as described at: http://www.nature.com/nbt/journal/v31/n9/full/nbt.2647.html

For a given nucleotide sequence, aggregate_score returns a number in (0,1] representing the risk for
off-target effects were that sequence to be used as a guideRNA for Cas9. 0 represents *high* risk
and 1 represents *low* risk. Therefore we can consider the aggregate_score as a measure of QUALITY.

aggregate_score is based on the summation of individual scores. So, if Blast search for a candidate guideRNA
returns 5 off-target sites, the aggregate_score will be 1/(1 + sum of individual off-target scores)

Individual scores are calculated by single_hit_score. For a pair of sequences, the single_hit_score combines
the total number of mismatches between those sequences, the positions of the mismatches, and the pairwise
distance between mismatches. A *low* single_hit_score corresponds to a low likelihood of off-target effects
between the two nucleotide sequences (i.e. low similarity), and therefore a *lower* aggregate score. Therefore
we can consider the single_hit_score as a measure of SIMILARITY, e.g. RISK.

The single_hit_score iterates through the mismatch positions between the two sequences. The aggregate score starts
at a maximum of 1. For a mismatch at position i, the score is multiplied by M[i], where M is an 
experimentally-determined array of mismatch penalties. The aggregate score is also multiplied by (1/n^2),
where n is the number of mismatches; this serves as a damping penalty for highly mismatched sequences.
Finally, the aggregate score is multiplied by the mean pairwise distance between mismatch sites. Thus two 
adjacent mismatches are considered less 'risky' than two distant mismatches.
"""

M = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]

def count_mismatch(guide,offtarget):
	assert len(guide) == len(offtarget)
	counter = 0
	for idx in range(len(M)):
		if guide[idx] != offtarget[idx]:
			counter+=1
	return counter


def mismatch_positions(guide,offtarget):
	assert len(guide) == len(offtarget)
	return [i for i in range(len(guide)) if guide[i]!=offtarget[i] and i<len(M)]


def single_hit_score(guide,offtarget):
	n_mismatch, pos_mismatch = count_mismatch(guide,offtarget), mismatch_positions(guide,offtarget)
	assert n_mismatch == len(pos_mismatch)
	pairwise_mean = mean_pairwise_distance(pos_mismatch)
	pairwise_damping = 1.0/(((19-pairwise_mean)/19.0)*4 + 1)
	try:
		mismatch_damping = 1.0/(n_mismatch**2)
	except ZeroDivisionError:
		mismatch_damping = 1
	final_score = 1
	for i in pos_mismatch:
		final_score *= (1-M[i])
	final_score*=(pairwise_damping*mismatch_damping*100)
	return final_score


def mean_pairwise_distance(iterator):
	if len(iterator) < 2:
		return 1
	pairwise_distances = []
	for i in range(len(iterator)):
		for j in range(i,len(iterator)):
			if i!=j:
				pairwise_distances.append(abs(iterator[i]-iterator[j]))
	return sum(pairwise_distances)/len(pairwise_distances)


def aggregate_score(guide, list_of_offtargets):
	total_score = 0
	if not list_of_offtargets:
		return 1
	for offtarget in list_of_offtargets:
		#print single_hit_score(guide, offtarget)
		total_score += single_hit_score(guide, offtarget)
	return 100.0/(100 + total_score)


if __name__ == "__main__":
	print len(M)
	test_guide = "GACGGGCGGCGGATCGGCAAAGG"
	offtargets = """GCCCGACGGGGGATCGGCAAGGG
GTGGCGCGGCGGATGGGCAAAGG
GATGGGCTGAGGATGGGCAATGG
GAGGGGCTGCAGTTCGGCAAGAG
GTCTGGCGGGGGATCGGGAAGGG
AACGAGCTGCGGATCTGCAATAG
GATGGGCTGCGGATGGGCACAGG
GCCAGGCGGCCGAGCGGCAACAG
GATGGGCGGCGGAGTGGCAAGGG
GGCGGGCGGCGGGTCGCCATCGG
GACGGGCTGGGGATCTGCAGCGG
GACGGGGAGCTGATCTGCAATGG
GAGGGGCGGGGGAGCGGCAGCGG
GGCGGGCGGCGGGTCGGCCGTGG
GTCGGGCGGGGGGGCGGCAAGAG
GAGGGGCGGGGGATAGGAAAAGG
GATGGGCGGCGGCTGGGCATTGG
GACTGGGGGCGGGTCCGCAAAAG
GAAGGGCGGAGGATTAGCAAGAG
GACAGGCGGCAGTTCAGCAATAG
GAGGGGCGGGGGAGAGGCAAAAG
GTCGGGCGGCGGCTGGGCCATGG
GGCGGGCGCCGGGGCGGCAAAGG
GACGGGCGGGCGAGCGGCAGGAG
GACGGGTGCCGGATCGGGTAGGG
GGCGGGCGGCGGTTTGGAAAAGG
GACGCGCGGCGGGGCGGCGAGAG
GACGGGGGGCGAACCGGCGAGGG
GCCGGGCGGCGGAAAGGAAACGG
GACGGGCGGTGGAACACCAAAGG
GACGGGCGGCAGATCCGAATGGG""".split()
	for off in offtargets:
		print off, mismatch_positions(off, test_guide), single_hit_score(test_guide, off)
	print aggregate_score(test_guide,offtargets)
