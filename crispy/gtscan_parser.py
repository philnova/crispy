gt = __import__('gtscan')
bowtiePath = '/Users/philnova/github_repos/crispy_site/bowtie-1.1.1/bowtie'
import os

def createCandidateSequenceString(candidate):
    candidateSequence = candidate.translate(None, '1234567890 ')
    rcCandidateSequence = gt.Seq(candidateSequence,gt.generic_dna).reverse_complement().tostring()
    return (candidateSequence, rcCandidateSequence, 'faName')

def findCandidateSequenceInfo(candidate,refGenome):
    candidateSequence,rcCandidateSequence,faName = createCandidateSequenceString(candidate)
    candidateSequenceLength = len(candidateSequence)
    if candidateSequenceLength<200:
        call= '{0} {1} -c {2}'.format(bowtiePath,refGenome,candidateSequence)
    else:
        call= '{0} {1} -c {2}'.format(bowtiePath,refGenome,candidateSequence[0:100])
    proc = gt.sp.Popen(call,shell=True,stderr=gt.sp.PIPE,stdout=gt.sp.PIPE)
    stdout_value, stderr_value= proc.communicate()
    zeroBasePosition = stdout_value.split('	')[1:4] # Results in array like ['+', '21', '33031934']
    try:
        candidateSequenceLocation = [zeroBasePosition[0],zeroBasePosition[1],int(zeroBasePosition[2])+1,candidateSequenceLength]
    except IndexError:
        candidateSequenceLocation = ['','','',candidateSequenceLength]
    return candidateSequenceLocation, faName

def findCandidateTargets(rule,candidate):
    print('Locating candidate targets...\n')
    targetRe = gt.createCandidateTargetRegex(rule)
    candidateSequence, rcCandidateSequence, faName = createCandidateSequenceString(candidate)
    windowSize = len(rule)
    candidateTargets = {}
    negFix = len(candidateSequence)-len(rule)


    # Scans string for candidate targets
    for i in range(0,len(candidateSequence)-windowSize+1):               ## Sliding window of len(rule)
        if targetRe.match(candidateSequence[i:i+windowSize]) != None:    ## If there IS match with regex 
            GC=gt.calculateGC(candidateSequence[i:i+windowSize])                 ## Also calculate GC content
            candidateTargets[i+1] = (targetRe.match(candidateSequence[i:i+windowSize]).group(),i+1,GC,'+')

    # Scans reverse complement string for candidate targets
    for i in range(0,len(rcCandidateSequence)-windowSize+1):               ## Sliding window of len(rule)
        if targetRe.match(rcCandidateSequence[i:i+windowSize]) != None:    ## If there IS match with regex 
            GC=gt.calculateGC(rcCandidateSequence[i:i+windowSize])                 ## Also calculate GC content
            candidateTargets[negFix-i+1+100000] = (targetRe.match(rcCandidateSequence[i:i+windowSize]).group(),negFix-i+1,GC,'-')
    return candidateTargets

#print createCandidateSequenceString('AGGCCCAAAAACCCTTAGG')

def find_offtargets(candidate, genome):
    
    offtargetfilter = "N" * (len(candidate) - 3) + "NRG"
    rule = "x" * (len(candidate) - 8) + "XXXXXNGG"
    mismatches = 2
    iupacdict = {
    'M':'[AC]',
    'R':'[AG]',
    'W':'[AT]',
    'S':'[CG]',
    'Y':'[CT]',
    'K':'[GT]',
    'V':'[ACG]',
    'H':'[ACT]',
    'D':'[AGT]',
    'B':'[CGT]',
    'N':'[ACGT]',
    'G':'G',
    'A':'A',
    'T':'T',
    'C':'C',
    }
    config = gt.ConfigParser.SafeConfigParser()
    config.readfp(open(gt.os.path.join(gt.base_path, 'config.ini')))
    refGenomeDir = config.get('important','ref_genome_dir')

    filterRegexString = ""
    for n in offtargetfilter:
        filterRegexString+=iupacdict[n]
    filterRegex = gt.re.compile(filterRegexString)
    
    
    # This is to parse the genome list from config.ini to work with Ensembl
    try:
        refGenomeInfo = config.get('ref_genomes',genome)
        ucscVars = refGenomeInfo.split('|')[1]
    except:
        ucscVars = genome

    
    refGenomeDir = config.get('important','ref_genome_dir')

    #important functions that do stuff
    whereabouts, faName = findCandidateSequenceInfo(
        candidate,gt.os.path.join(refGenomeDir,genome))
    candidateTargets = findCandidateTargets(
        rule, candidate)
    potentialOffTargetStrings = gt.createPotentialOffTargetStrings(
        rule, candidateTargets)
    offTargSummary,potentialOffTargets = gt.filterPotentialOffTargets(
        candidateTargets, potentialOffTargetStrings, os.path.join(refGenomeDir,genome), mismatches, rule, filterRegex, mismatches, whereabouts)
    return potentialOffTargets

if __name__ == "__main__":
	print find_offtargets('GACGGGCGGCGGATCGGCAAAGG','GCA_000001405.15_GRCh38_no_alt_analysis_set')



    