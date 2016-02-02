import sequence as seq
import sys
from models import Sequence

def count_lines(filename):
	with open(filename) as fo:
		for idx, line in enumerate(fo):
			pass
	return idx

def output(percent_complete):
	sys.stdout.write('%s\r' % percent_complete)
	sys.stdout.flush()

def offtarget_score(sequence):
	return .5

def main(chrm_name, limit=float('inf')):
	filename = 'MergedGuideRNA/chr'+str(chrm_name)+'.txt'
	count = count_lines(filename)
	print "{0} lines in file".format(count)
	with open(filename) as fo:
		for idx, line in enumerate(fo):
			if not idx%1000:
				output("{0:.2f}".format(float(idx) * 100/count) + ' percent complete')
			if idx<limit and idx:
				new = seq.NewSequence()
				linelist = line.split()
				chromosome, start, stop, sequence = str(chrm_name), int(linelist[1]), int(linelist[2]), linelist[3]
				pam = sequence[-3:]
				score = offtarget_score(sequence)
				if linelist[-1] == 'F':
					direction = 'forward'
				else:
					direction = 'reverse'
				new.chromosome, new.start, new.stop, new.sequence, new.direction, new.organism, new.assembly, new.pam, new.score = chromosome, start, stop, sequence, direction, 'human', 'hg19', pam, score
				new.load()

def all_chromosomes():
	for i in xrange(1,23):
		main(i)

if __name__ == "__main__":
	main('Y',1000)
	for i in seq.session.query(Sequence).limit(10):
		print i.sequence, i.start, i.stop