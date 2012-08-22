#!/usr/bin/env python3
# codonopt.py - A codon-optimising script with support for enzyme site exclusion.
# Uses a weighted-random codon selection system with optional minimum frequency threshold to exclude very rare codons.
# Comes with presets including empirically determined E.coli optimal codon tables from Welch et al, 2009.
# Applies miscellaneous rules from "Designing Genes for Successful Protein Expression" by Welch et al, 2011:
#  - TODO: Minimum frequency threshold is increased for first 25 codons.
#  - TODO: "NGG" codons are penalised in first 25 codons.

import random
import sys
import os
import argparse
import json
from math import ceil, floor

ScriptDir = sys.path[0] # Gives script directory. Default/imported codon tables will be here.
LibDir = os.path.join(ScriptDir, 'lib')
CodonTableDir = os.path.join(LibDir, 'codontables')
WorkingDir = os.getcwd() # Gives current working directory.

ArgParser = argparse.ArgumentParser(description=('Performs reverse-translation and weighted-random '
                                         'codon optimisation on an amino acid sequence. '
                                         'Default setting is an E.coli table derived from '
                                         'empirical synthetic DNA research by Welch et al, '
                                         '2009.'))
ArgParser.add_argument('infile', type=str,
                       help='Protein sequence, in FASTA file format, to be reverse translated and codon optimised.')
ArgParser.add_argument('-d', '--output-dna', action='store_true',
                       help='Output codon-optimised sequence as DNA, not RNA.')
ArgParser.add_argument('-v', '--verbose', action='store_true',
                       help='Print far more information during processing.')
ArgParser.add_argument('-s', '--species', type=str, default='default_ecoli',
                       help='Species to optimise DNA for. Default is E.coli.')
ArgParser.add_argument('-x', '--exclude-sites',
                       help=('A json-formatted list of DNA or RNA substrings to forbid, '
                             'usually restriction enzyme target sites.'))
ArgParser.add_argument('-m', '--min-codon-frequency', type=float, default=0.0,
                       help=('A minimum codon frequency to accept; frequencies below this are discarded and'
                            ' remaining codon frequencies re-normalised.'))
ArgParser.add_argument('-n', '--splice-candidates', type=int, default=0,
                       help=('Number of random candidate codon sets to generate before attempting to splice '
                             'together a suitably enzyme-free set. Larger numbers of candidates may be useful '
                             'when attempting to exclude many sites, or highly redundant sites.'))
ArgParser.add_argument('--usage', action='store_true',
                       help='Ignore other flags/input and provide detailed usage information.')
# Future features:
#ArgParser.add_argument('--convert-table', type=FileType('r'),
#                       help='Target table in XYZ format for import to native JSON format')

Args = vars(ArgParser.parse_args()) # vars makes a dict of the results.
# ============================================================================================== #
def verboseMsg(Msg): # Optional print statement.
    if Args['verbose']: print(Msg)

class RevTranslateRNA:
    def __init__(self, aminos, codontable, excludesites, MinCodonFrequency=0.0, outputdna=False, candidates=0):
        '''aminos: A string containing just the peptide sequence of a desired protein.
           codontable: Should be a dict mapping Amino acid letters to codons with frequencies. Example entry:
              "C": { "TGC": { "frequency": 0.5772005772005773 }, 
                     "TGT": { "frequency": 0.42279942279942284}  }
           exclude_sites: Should be a list of sites to exclude; list is iterated over and compared to output.
           output_dna: Boolean: output is converted to DNA before export.
           candidates: Number of permutations of raw codon list to create by weighted-random selection before splicing.
        '''
        def ErrThenExit(String):
            print(String)
            sys.exit(1)

        try: self.peptidesequence = self.loadaminos(aminos)
        except AssertionError: ErrThenExit("Error: Provided 'aminos' not a string!")

        try: self.excludesites, self.largestexclude = self.parseExcludes(excludesites)
        except AssertionError: ErrThenExit("Error: Provided 'excludesites' not a list of strings!")

        try: self.codontable = self.importCodonTable(codontable, MinCodonFrequency)
        except AssertionError: ErrThenExit("Error: Provided 'codontable' not a dictionary!")
        self.translationtable = self.buildTranslationTable(codontable)

        self.buildwindowsize = ceil(self.largestexclude/3) # i.e. round up.
        self.outputdna = outputdna
        self.rnacodonlist = []
        self.ignoredCodons = [] # Used to ignore intransigent sites.
        self.result = ''
        candidatesnum = candidates or 20
        self.generateCandidateList(candidates)
        verboseMsg("Initialised RevTranslateRNA object.")

    def translate(self,codonlist):
        'Translate a list of codons into a string of amino acids. No regard for start/stop: translates all.'
        outstring = ''
        for codon in codonlist:
            encodes = self.translationtable[codon.upper()]
            outstring += encodes
        return outstring

    def reverseTranslate(self, AminoString):
        'Returns a list of codons encoding the AminoString, chosen by weighted-random from self.codontable.'
        CodonList = []
        for Amino in AminoString:
            Codon = self.selectCodon(Amino)
            CodonList.append(Codon)
        return CodonList

    def generateCodonList(self):
        return self.reverseTranslate(self.peptidesequence)

    def generateCodonSequence(self):
        'Shortcut method that just sets self.rnacodonlist to a raw reverse-translated form of self.peptidesequence.'
        self.rnacodonlist = self.generateCodonList()

    def generateCandidateList(self, nums=20):
        assert nums > 0 # Better an assertionerror than an infinite loop on negative nums..
        self.candidatecodonlists = []
        while nums > 0:
            nums -= 1 # Count down.
            self.candidatecodonlists.append(self.generateCodonList())

    def prettyCodonsList(self, Codons):
        'An improvement to the codon output format; tabulates codons in rows of ten for ease of manual error checking.'
        NumCodons = len(Codons)
        FormatString = '{}\t{} {} {} {} {} {} {} {} {} {}\r\n'
        PrintString = ''
        for x in range(0, NumCodons, 10): # x counts from 0 in increments of 10
            sublist = Codons[x:x+10]
            while len(sublist) < 10: # Prevent format errors on last line.
                sublist.append('')
            PrintString += FormatString.format(x,*sublist)
        return PrintString

    def spliceSequences(self, SequenceList=[]):
        '''Takes a list of candidate sequences (which must encode same amino sequence) and attempts to splice a site-free form.
        
        This method uses self.issueMap to do the job of mapping each candidate sequence for issues, and compares sequences
        to find which one has the longest stretch of codons without any issues to address. It appends all the codons prior to
        the first discovered issue to an internal candidate sequence, and then starts again from this point, attempting to find
        another candidate with a useful stretch of site-free nucleotides.

        In doing so, two problems may arise:
            A) The method may encounter an issue/site that is universal to all sequences, and it must ignore said site and continue.
            B) The method may create new issues/sites at the splice junctions, but resolving these issues can be left to the
               self.optimiseSequence() method.
        '''
        SeqLen = 0 # Changed by assertEquivalence to length of sequences.
        def assertEquivalence(SeqList):
            'This function checks that SeqList is a non-zero list and that contents are valid and equivalent codon lists.'
            assert isinstance(SeqList, list)
            assert len(SeqList) > 0
            for i in SeqList:
                assert isinstance(i, list)
            First = self.translate(SeqList[0])
            nonlocal SeqLen
            SeqLen = len(SeqList[0]) # Used later for a while loop.
            for i in SeqList:
                assert self.translate(i) == First
            return True

        if not SequenceList:
            SequenceList = self.candidatecodonlists
        try: assertEquivalence(SequenceList)
        except AssertionError:
            verboseMsg("Splicer: ERROR: SequenceList must be a list of codon lists, which all encode the same amino sequence.")
            verboseMsg("Splicer: ERROR: Can't work under these conditions, returning None.")
            return None
        
        # Below: We're going to iterate over SequenceList, making a new working list with dictionary values instead of plain codon lists.
        # Each dictionary is simply the issueMap of the codon list, with the codon list embedded into it under key "fullcodonlist".
        WorkingList = []
        for i in SequenceList:
            i_issueMap = {}
            i_issueMap['issueMap'] = self.mapExcludedSites(i, self.excludesites)
            i_issueMap['fullcodonlist'] = i
            i_issueMap['mapIndices'] = list(i_issueMap['issueMap'].keys())
            i_issueMap['mapIndices'].sort() # Probably unnecessary, as dict-keys are already sorted..
            WorkingList.append(i_issueMap)
            # i_issueMap['issueMap'] is a dict with numerical keys.
            # For each key (which corresponds to a site with issue(s)), the values returned are:
            #   {'spanning_codons', 'encodes', 'contains_sites', 'ensuing_codons', 'preceding_codons'}
            # Splicer will (probably) ignore most of these keys, instead comparing the 'mapIndices' lists simply to check which
            # Sequence has the longest span of uninterrupted codons.

        NewCodons = []
        verboseMsg("Splicer: About to attempt splicing of the following sequences:")
        for i in WorkingList:
            i_number = WorkingList.index(i)
            verboseMsg("Splicer: Seq #"+str(i_number)+" issues: "+str(i['mapIndices']))

        def findBestSpan(WorkList, CodonOffset):
            'Working from CodonOffset, compare each sequence in WorkList to find one with longest-stretch-til-next-issue.'
            verboseMsg("Splicer: Working from Offset "+str(CodonOffset)+"...")
            BestSoFar = {'Seq':-1, 'Codon':-1} # Index of which sequence has best stretch of clear codons.
            for i in range(0,len(WorkList)):
                i_indices = WorkList[i]['mapIndices']
                for x in i_indices[:]: #Scan over a slice-copy of the list so we can edit the list as we iterate.
                    if x < CodonOffset: del(i_indices[i_indices.index(x)]) # Delete if X less than CodonOffset.
                if len(i_indices) < 1:
                    BestSoFar = {'Seq':i, 'Codon':SeqLen-1}
                    continue
                if i_indices[0] > BestSoFar['Codon']:
                    BestSoFar = {'Seq':i, 'Codon':i_indices[0]}
            verboseMsg("Splicer: Best so far is Seq #"+str(BestSoFar['Seq'])+", which spans up until codon "+str(BestSoFar['Codon'])+".")
            if (BestSoFar['Codon'] == CodonOffset) or (BestSoFar['Seq'] < 0):
                verboseMsg("Splicer: Hit a brick wall at codon #"+str(CodonOffset)+" (recalcitrant enzyme site, or last codon), skipping.")
                return WorkList[BestSoFar['Seq']]['fullcodonlist'][CodonOffset:CodonOffset+1] #Return next codon of "best" anyway.
            BestFullLength = WorkList[BestSoFar['Seq']]['fullcodonlist']
            NextCodons = BestFullLength[CodonOffset:BestSoFar['Codon']]
            #verboseMsg("Splicer: Returning codons for extension:\n"+str(NextCodons))
            return NextCodons

        while len(NewCodons) < SeqLen:
            NewCodons += findBestSpan(WorkingList, len(NewCodons))
        verboseMsg("Splicer: Sanity checking that NewCodons matches amino sequence of input sequence set..")
        assert self.translate(NewCodons) == self.translate(SequenceList[0])
        self.rnacodonlist = NewCodons

    def optimiseSequence(self):
        import itertools
        def createCodonPermutations(AminoSequence, sanity_ceiling = 256):
            '''Return a list of codon-lists that encode the desired AminoSequence. Warning: for long
            sequences, this might be a really dumb thing to do, so this implements a sanity check,
            comparing the predicted number of permutations to the "sanity_ceiling" and aborting if
            greater.'''
            CodonChoices = [] # Will be populated with sub-lists of codon choices at every index.
            for Amino in AminoSequence:
                PossibleCodons = list(self.codontable[Amino.upper()].keys())
                CodonChoices.append(PossibleCodons)
            PredictedPermutations = 1
            for CodonSet in CodonChoices: # Calculate # of possible permutations by multiplying choices
                PredictedPermutations *= len(CodonSet)
            if PredictedPermutations > sanity_ceiling:
                verboseMsg("Predicted permutations = "+str(PredictedPermutations)+".. aborting!")
                return []
            verboseMsg("Reverse Translator: Generating "+str(PredictedPermutations)+" permutations for amino sequence "+AminoSequence+"..")
            PermutationMaker = itertools.product(*CodonChoices)
            Perms = []
            for Permutation in PermutationMaker:
                Perms.append(Permutation)
            return Perms
        def chooseBestCodonPermutation(Permutations):
            'This is here to be improved upon.'
            import random
            return random.choice(Permutations)
        def shuffleCodons(CodonList):
            'Translates and then reverse-translates the target region.'
            Aminos = self.translate(CodonList)
            NewCodons = self.reverseTranslate(Aminos)
            return NewCodons
        def attemptPermutations(IssueDict):
            IssueContext = ' '.join(IssueDict['preceding_codons']) + '{0}' + ' '.join(IssueDict['ensuing_codons'])
            Permutations = createCodonPermutations(IssueDict['encodes'])
            if Permutations: # i.e. if we were not handed back an empty list..
                HitList = [] # List permutations to remove.
                verboseMsg("Reverse Translator: Trimming permutations to remove those containing excluded sites..")
                verboseMsg("Reverse Translator: Permutations are: "+json.dumps(Permutations))
                for Perm in Permutations:
                    PermSequence = IssueContext.format(''.join(Perm)).upper() # Insert Permutation into local context..
                    for Exclude in self.excludesites: # Then scan through excludes..
                        if Exclude in PermSequence: # ..and remove the permutation if an excluded site is found within.
                            HitInd = Permutations.index(Perm)
                            HitList.append(Perm)
                            verboseMsg("Reverse Translator: Targeting permutation index '"+str(HitInd)+"' for deletion.")
                for Hit in HitList:
                    try: Permutations.remove(Hit)
                    except ValueError: verboseMsg("Permutation "+str(Hit)+" already removed from permutations.")
                verboseMsg("Reverse Translator: Permutations are: "+json.dumps(Permutations))
                if len(Permutations) < 1: # if there are no options left!
                    print("Reverse Translator: Could not find an enzyme-free permutation for issue at index "+str(Index)+", "+\
                          "corresponding to nucleotide "+str(Index*3)+"! You must address this manually. Sorry!")
                    self.ignoredCodons.append(Index)
                else: # Randomly select a suitable codon set, and replace the problem area with it:
                    SuitableCodons = chooseBestCodonPermutation(Permutations) #Currently random.
                    NumCodons = len(SuitableCodons)
                    verboseMsg("Reverse Translator: Choosing: "+''.join(SuitableCodons) + " for codons "+str(Index)+" to "+str(Index+NumCodons-1)+".")
                    verboseMsg("New Context is:\r\n"+IssueContext.format('|'+' '.join(SuitableCodons)+'|').upper())
                    self.rnacodonlist[Index:Index+NumCodons] = SuitableCodons
            else: # If we were given an empty list, probably because number of permutations was too large.
                print("Reverse Translator: There are unresolved issues at index "+str(Index)+", corresponding to nucleotide "+\
                       str(Index*3) + ". Due to the number of permutations at this locus, it has been skipped."+\
                       " You must manually edit this site to your satisfaction. Sorry!")
                self.ignoredCodons.append(Index)

        verboseMsg("Reverse Translator: Generating first-pass codon list for amino acid sequence.")
        if not self.rnacodonlist: # If there's nothing to work with, generate a candidate:
            self.rnacodonlist = self.reverseTranslate(self.peptidesequence)
            verboseMsg("Reverse Translator: First pass codon list generation complete.")
        i = 3 # Really it should be done by 3 cycles..
        while i > 0:
            i -= 1
            issueMap = self.mapExcludedSites(self.rnacodonlist, self.excludesites)
            if not issueMap: #i.e. if no issues remain
                break

            # issueMap is a dict with numerical keys. The values returned are:
            #   {'spanning_codons', 'encodes', 'contains_sites', 'ensuing_codons', 'preceding_codons'}
            # Procedure: 1) Generate all permitted permutations of codons that encode "encodes" (including preceding amino)
            #            2) Scan permutations in their proposed context and delete any that contain anything in "contains_sites".
            #            3) From remainder, choose randomly.

            IndexList = list(issueMap.keys()) # Create an ordered list of indices. (for some reason, using sort() yielded Nonetype?
            verboseMsg("Reverse Translator: IndexList is: "+str(IndexList))
            for Index in IndexList:
                verboseMsg("Reverse Translator: Addressing issue at codon "+str(Index)+" (nucleotide "+str(Index*3)+"):")
                ProblemDict = issueMap[Index]
                attemptPermutations(ProblemDict)

        # Finally:
        verboseMsg("Reverse Translator: Finished addressing issues, concatenating codons to result!")
        self.result = ''.join(self.rnacodonlist)
        if self.outputdna:
            self.result = self.result.replace("U", "T")

    def mapExcludedSites(self, CodonList, ExcludedSites):
        '''Returns a dictionary containing integer keys, where each key corresponds to a list index for
        CodonList where an excluded site was found. The dictionary contains the number of codons the
        excluded site spans, the amino acids encoded by those codons, and any other sub-sites found at
        the same index.'''
        def allindices(string, sub, offset=0):
            'Returns a list of ALL indices of a substring, not just first.'
            listindex = []
            i = string.find(sub, offset)
            while i >= 0: # Because find returns -1 when it finds nothing.
                listindex.append(i)
                i = string.find(sub, i + 1)
            return listindex
        def findProblemCodons(sequence, excludedsite):
            'Maps string indices to codon-list indices.'
            stringindices = allindices(sequence, excludedsite)
            codonlistindices = []
            for stringindex in stringindices:
                codonindex = floor(stringindex/3) # floor = imported from math.
                if codonindex not in self.ignoredCodons:
                    if codonindex not in codonlistindices:
                        codonlistindices.append(codonindex)
                else:
                    verboseMsg("Mapper: Ignoring issue at index "+str(codonindex)+" because it can't be resolved.")
            return codonlistindices

        verboseMsg("Mapping issues in RNA codon map.")
        verboseMsg(self.prettyCodonsList(CodonList))

        FoundIssues = {}
        sequenceString = ''.join(CodonList)
        for Exclude in ExcludedSites:
            ExcludeCodonLength = ceil(len(Exclude) / 3) # Round up.
            ProblemCodons = findProblemCodons(sequenceString, Exclude)
            for ProblemCodon in ProblemCodons: # ProblemCodon is a list index for the codon in question.
                verboseMsg("Mapper: Found site: "+Exclude+" @ codon "+str(ProblemCodon)+".")

                # Below: Sites don't always align perfectly with boundary of "problemcodon", so must determine if overhang exists.
                #  So, join codons, search for problem site, and if it's not in joined codons, probably need extra codon.
                if Exclude in ''.join(CodonList[ProblemCodon : ProblemCodon + ExcludeCodonLength]):
                    SpannedCodons = CodonList[ProblemCodon : ProblemCodon + ExcludeCodonLength]
                    Ensuing = ProblemCodon + ExcludeCodonLength
                else:
                    SpannedCodons = CodonList[ProblemCodon : ProblemCodon + ExcludeCodonLength + 1]
                    Ensuing = ProblemCodon + ExcludeCodonLength + 1
                    assert Exclude in ''.join(SpannedCodons) # Throw a fit if this fails.

                PrecedingCodons = CodonList[ProblemCodon-self.largestexclude:ProblemCodon]
                EnsuingCodons = CodonList[Ensuing:Ensuing+self.largestexclude]

                ProblemString = '\tContext: {0}|{1}|{2}'.format(' '.join(PrecedingCodons), '-'.join(SpannedCodons), ' '.join(EnsuingCodons))
                verboseMsg(ProblemString)

                if ProblemCodon not in FoundIssues.keys():
                    FoundIssues[ProblemCodon] = {'spanning_codons':	ExcludeCodonLength,
                                                'encodes': 		self.translate(SpannedCodons),
                                                'preceding_codons':	PrecedingCodons,
                                                'ensuing_codons':	EnsuingCodons,
                                                'contains_sites': 	[Exclude]              }
                elif FoundIssues[ProblemCodon]['spanning_codons'] >= ExcludeCodonLength:
                    FoundIssues[ProblemCodon]['contains_sites'].append(Exclude)
                elif FoundIssues[ProblemCodon]['spanning_codons'] < ExcludeCodonLength:
                    FoundIssues[ProblemCodon]['spanning_codons'] = ExcludeCodonLength
                    FoundIssues[ProblemCodon]['contains_sites'].append(Exclude)
                    FoundIssues[ProblemCodon]['encodes'] = self.translate(SpannedCodons)
                else:
                    print("Error in mapExcludedSites.")
        verboseMsg("Mapper: Finished! Returning to reverse translator.")
        verboseMsg("Mapper: Issues are: "+json.dumps(FoundIssues,indent=1))
        return FoundIssues

    def selectCodon(self, Amino, Floor=0.0):
        '''Returns a weighted-random-selected codon with frequency greater than "floor".'''
        def randomCategory(probDict):
            '''Accepts a dictionary of {opt:wgt}, e.g.: randomCategory({'a':.15, 'b':.35, 'c':.5})
               Returns a selection based on weight. Weights must be normalised, should sum to 1!'''
            import random
            r = random.random() # range: 0,1
            total = 0
            for value,prob in probDict.items():
                if prob <= 0: continue # ignore items with a 0 weight.
                total += prob
                if total>r: return value
            raise Exception('Distribution not normalized: {probs}'.format(probs=str(probDict)))
        # Above function is re-usable. Below wraps it for specific codon-selection purposes.
        assert isinstance(Amino, str)
        CandidateCodons = self.codontable[Amino]
        ChosenCodon = randomCategory(CandidateCodons)
        return ChosenCodon

    def parseExcludes(self, ExcludeList):
        verboseMsg("Parsing through excluded enzyme list.")
        def cleanupSite(Site):
            CharHitList = '0123456789()<>^\\/-[],.'
            for Char in CharHitList:
                Site = Site.replace(Char, '')
            Site = Site.strip().strip("N") # Kill whitespace and random flanking Ns!
            return Site.upper()
        assert isinstance(ExcludeList, list)
        OutputList = []
        Largest = 0
        for entry in ExcludeList:
            assert isinstance(entry, str)
            Entry = cleanupSite(entry)
            EntryLen = len(Entry)
            OutputList.extend(self.createRNAPermutations(Entry))
            if EntryLen > Largest:
                Largest = EntryLen
        verboseMsg("Excluded sites are: \r\n"+str(OutputList))
        return OutputList, Largest            

    def buildTranslationTable(self, codontable):
        'Inverts a codon table, omitting frequencies, to infer the translation table of the target species.'
        TransTable = {}
        verboseMsg("Building forward translation table from codon table.")
        for Amino in codontable.keys():
            for Codon in codontable[Amino]:
                TransTable[Codon] = Amino
        return TransTable

    def loadaminos(self, aminos):
        assert isinstance(aminos, str)
        # NOTE: This function merely strips illegal (non-IUPAC) characters.
        #  It expects JUST the amino sequence, not markup e.g. FASTA titles! Sanitise input!
        permittedaminos = 'ACDEFGHIKLMNPQRSTVWY' # TODO: implement B, X and Z?
        aminos = aminos.upper()
        for char in aminos:
            if char not in permittedaminos:
                print("Found illegal character (not an IUPAC amino acid): "+char+"..deleting all instances.")
                aminos = aminos.replace(char, "")
        return aminos

    def importCodonTable(self,CodonTable, Floor = 0.0):
        '''Imports codon tables that may have excess metadata into simple codon:frequency dicts.
           Also removes codons whose frequency is below Floor, and re-normalises remaining frequencies.
           i.e. convert this:
              "C": { "TGC": { "frequency": 0.5772005772005773, 'info':'More highly charged under starvation' }, 
                     "TGT": { "frequency": 0.42279942279942284, 'synthetase':{'seq':'mabcd','protID=foo1234'}
                    }
           ...into this:
              "C":{"TGC":0.5772005772005773,
                   "TGT":0.42279942279942284}'''
        verboseMsg("Importing specified codon table, discarding frequencies below "+str(Floor))
        assert isinstance(CodonTable, dict)
        FrequencyDict = {}
        for Amino in CodonTable.keys():
            FrequencyDict[Amino] = {}
            ThisAminoTotalFrequency = 0.0
            for Codon in CodonTable[Amino].keys():
                # This For-loop gets the frequency for each codon, and copies it if it's greater than Floor.
                try:
                    Freq = CodonTable[Amino][Codon]['frequency']
                    if Freq > Floor:
                        FrequencyDict[Amino][Codon] = Freq
                        ThisAminoTotalFrequency += Freq
                except:
                    print("Error while importing codon table; entry {0} in {1} \
                           has no frequency value. Skipping {0}.".format(Codon,Amino))
            if ThisAminoTotalFrequency < 1:
                # This probably means we've omitted one or more codons due to being less frequent than floor.
                # So, we have to re-normalise the remaining frequencies:
                for Codon in FrequencyDict[Amino].keys():
                    FrequencyDict[Amino][Codon] = ThisAminoTotalFrequency / FrequencyDict[Amino][Codon]
        return FrequencyDict

    def createDNAPermutations(self, Site):
        'Returns a list of possible permutations for a given IUPAC DNA string.'
        verboseMsg("Creating a set of DNA permutations for "+Site)
        import itertools
        IUPACDNA = { 'A': 'A',    'B': 'CGT',
                     'C': 'C',    'D': 'AGT',
                     'G': 'G',    'H': 'ACT',
                     'K': 'GT',   'M': 'AC',
                     'N': 'ACGT', 'R': 'AG',
                     'S': 'GC',   'T': 'T',
                     'V': 'ACG',  'W': 'AT',
                     'Y': 'CT',   '.': 'ACGT',
                     '-': 'ACGT' }
        SNPList = []
        PermList = []
        for Char in Site:
            SNPList.append(IUPACDNA[Char])
        PermGenerator = itertools.product(*SNPList) # "*" character dumps list contents as args.
        for item in PermGenerator:
            PermList.append(''.join(item))
        return PermList

    def createRNAPermutations(self, Site):
        'Returns a list of possible RNA permutations for a given IUPAC DNA string.'
        DNASites = self.createDNAPermutations(Site)
        verboseMsg("Converting DNA permutations to RNA.")
        RNASites = []
        for DNASite in DNASites:
            RNASites.append(DNASite.replace("T", "U"))
        return RNASites

    def convertRNAtoDNA(self,RNAinput):
        assert isinstance(RNAinput, str)
        DNAoutput = RNAinput.upper().replace("T","U")
        return DNAoutput

    def giveOutput(self):
        'Runs quick sanity checks and returns result.'
        newaminosequence = self.translate(self.rnacodonlist)
        #print(newaminosequence)
        if not self.peptidesequence == newaminosequence:
            print(("ERROR: Input peptide sequence and translated peptide sequence of output RNA/DNA do not match!\n")
                  ("Please report this to the developer on Github: https://github.com/cathalgarvey\n")
                  ("--Input--\r\n{0}\n".format(self.peptidesequence))
                  ("--Output--\r\n{0}".format(newaminosequence)))
            sys.exit(1)
        if self.outputdna:
            self.result += 'TAA'
        else:
            self.result += 'UAA'
        return self.result

# ======================================================================== #
# Functions:
from fasta_utils import ReadFastaFile

def importCodonTable(CodonTableFile):
    '''Expects a file containing a codontable in json format.

    codontable: Should be a dict mapping Amino acid letters to codons with frequencies. Example entry:
    "C": { "TGC": { "frequency": 0.5772005772005773 }, 
           "TGT": { "frequency": 0.42279942279942284}  }
    The file can actually contain other metadata for each codon, but should not contain keys other
    than codons for each amino acid; keys are assumed to correspond to codons only. Likewise, keys
    in the "master" dictionary are assumed to be amino acids and will be imported as such, although
    there is less risk of errors here as keys will only become relevant if called by a function/method.'''
    import json
    CodonTableFileDir = os.path.join(CodonTableDir, CodonTableFile + '.json')
    with open(CodonTableFileDir) as TableFile:
        CodonTable = json.loads(TableFile.read())
    assert isinstance(CodonTable, dict)
    ReturnTable = {}
    for Amino in CodonTable.keys(): # Assumes all keys are amino acids.
        ReturnTable[Amino] = {}
        for Codon in CodonTable[Amino].keys(): # Assumes all keys are codons.
            # What does this achieve, exactly? It just prunes codontable and verifies that each
            #  codon has a frequency value..useful?
            try: ReturnTable[Amino][Codon] = {'frequency':CodonTable[Amino][Codon]['frequency']}
            except: print("Codon "+Codon+" of Amino "+Amino+" has no 'frequency' key! Skipped.")
    return ReturnTable

def importExclusions(ExcludeFile):
    '''Expects a JSON formatted file containing a simple list of IUPAC DNA sequences to avoid.

    Sites are expanded into permutations and stored in the object for comparison to draft sequences.
    Providing too many sites with ambiguous characters from the IUPAC set may result in huge exclusion
    lists, potential program slowdown or crashing. For the love of all that's holy, strip unnecessary "N"
    characters from TypeIIs or Homing Endonucleases, where sites may often be written with lots of Ns outside
    the actual recognition area; this program will assume these are relevant and will generate 4^(number of Ns)
    permutations.'''
    import json
    if not isinstance(ExcludeFile, str):
        return []
    with open(ExcludeFile) as ExcludesFile:
        ExcludesList = json.loads(ExcludesFile.read())
    assert isinstance(ExcludesList, list)
    IUPACDNA = 'ABCDGHKMNRSTVWY.-'
    for Exclude in ExcludesList:
        for char in Exclude:
            assert char in IUPACDNA
    return ExcludesList

# ======================================================================== #
# Args contains: 'infile', 'output_dna', 'species', 'exclude_sites', 'min_codon_frequency', 'usage', 'splice_candidates'
# Usage: RevTranslatedRNA(aminos, codontable, excludesites, output_dna=False):
FastaDict = ReadFastaFile(Args['infile'])
AminoSequence = FastaDict['sequence']
if not AminoSequence:
    print("No valid Amino Acid sequence found in target file. Quitting.")
    sys.exit(1)
try: SpeciesCodonTable = importCodonTable(Args['species'])
except AssertionError:
    print("Codon table file specified is not a dict. Quitting.")
    sys.exit(1)
except KeyError:
    print("Encountered a missing dictionary key (Amino? Codon?) while parsing and pruning Codon Table. Quitting.")
    sys.exit(1)
try: SiteExclusionList = importExclusions(Args['exclude_sites'])
except AssertionError:
    print(("Either the specified file doesn't encode a list,"
           "or the list entries aren't all valid IUPAC DNA strings.\r\n"
           "NB: Permitted IUPAC DNA codes are 'ABCDGHKMNRSTVWY.-'"
           "\r\nProgram will quit."))
    sys.exit(1)
if Args['output_dna']: DNAArg = True
else: DNAArg = False

if Args['min_codon_frequency']: FreqArg = Args['min_codon_frequency']
else: FreqArg = 0.0

if Args['splice_candidates']: SplicesArg = Args['splice_candidates']
else: SplicesArg = 0

ReverseTranslator = RevTranslateRNA(AminoSequence, SpeciesCodonTable, SiteExclusionList, FreqArg, DNAArg, SplicesArg)
ReverseTranslator.spliceSequences()
ReverseTranslator.optimiseSequence()

# Deliver results of algorithmic goodness:
verboseMsg("Codons with remaining issues: "+json.dumps(ReverseTranslator.ignoredCodons))
verboseMsg("\r\nOutput:\r\n")
print(ReverseTranslator.giveOutput())

