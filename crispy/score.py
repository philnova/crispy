M = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]

def count_mismatch(guide,offtarget):
	assert len(guide) == len(offtarget)
	counter = 0
	for idx in range(len(guide)):
		if guide[idx] != offtarget[idx]:
			counter+=1
	return counter


def mismatch_positions(guide,offtarget):
	assert len(guide) == len(offtarget)
	return [i for i in range(len(guide)) if guide[i]!=offtarget[i]]

def single_hit_score(guide,offtarget):
	n_mismatch, pos_mismatch = count_mismatch(guide,offtarget), mismatch_positions(guide,offtarget)
	assert n_mismatch == len(pos_mismatch)
	pairwise_mean = mean_pairwise_distance(pos_mismatch)
	pairwise_damping = 1.0/(((19-pairwise_mean)/19.0)*4 + 1)
	mismatch_damping = 1.0/(n_mismatch**2)
	final_score = 1
	for i in pos_mismatch:
		final_score *= (1-M[i])*pairwise_damping*mismatch_damping
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
		print single_hit_score(guide, offtarget)
		total_score += single_hit_score(guide, offtarget) * 100
	return 100.0/(100 + total_score)


if __name__ == "__main__":
	pass
