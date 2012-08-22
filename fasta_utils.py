#!/usr/bin/env python3
'''Some functions for working with FASTA sequences as native types.'''

def ReadFastaFile(File):
    'Imports the contents of a file as a string, passes to ReadFasta() and returns output.'
    with open(File) as FastaFile:
        FastaDict = ReadFasta(FastaFile.read())
    return FastaDict

def ReadFasta(FastaString):
    '''Finds a sequence and returns a dictionary containing title, comments, accession (if found) and sequence.

    Supports extraction of titles, comments and sequence. (Comments are lines beginning with ';'.
    This is an oft-ignored FASTA spec feature and is almost never used.)
    Dictionary returned contains following keys, whether empty (evals false) or with contents:
     'title': (a string containing the title, including the ">" symbol)
     'comments': (a dict of format {string index: comment})
     'sequence': (a string containing the sequence)
     'accession': (a string containing the accession information if present (detected by searching for '|' in title)
    '''
    Sequence = {'title':'','comments':{},'sequence':'','accession':''}
    FastaLines = []
    for Line in FastaString.split("\n"):
        Line = Line.strip()
        if Line.startswith(">"): # Is sequence descriptor
            Sequence['title'] = Line[1:].strip()
        elif Line.startswith(";"): # Is a comment, rarely used
            Sequence['comments'][len(Sequence['sequence'])] = Line[1:].strip()
            # Uses sequence length SO FAR minus 1, so acts like a string index for comment markup location.
            # Sequence['comments'].append(Comment)
        else:
            Sequence['sequence'] += Line.upper().strip()
    FirstTitleWord = Sequence['title'].split(' ')[0].strip()
    if '|' in FirstTitleWord: # Yes, this is crude, I know. #lazy
        Sequence['accession'] = FirstTitleWord
        Sequence['title'] = Sequence['title'][len(Sequence['accession']):].strip()
    return Sequence

def MakeFasta(PassedSequenceDict, wraplines = 50):
    '''Returns a multiline string containing the FASTA equivalent of the passed dictionary.
    
    Dictionary must contain at least 'title' and 'sequence' keys. It may also contain 'accession'
    and 'comments' keys.

    This function throws an AssertionError if "sequence" or "title" keys are missing, or if not a
     dict, or if not all "comments" keys are integers.

    Accession should be a string, and will be prepended to output title.
    Comments should be a dictionary, with each key being the numeric string index where the
     comment is to occur, with the value being the comment, like so: {12:"Comment at position 12"}

    Comments are prepended with ';' as specced for FASTA, but are not compatible with many/most
    software packages.

    wraplines is the number of characters to print per line before wrapping to next line.
    spacelines is optional: place a space every (spacelines) characters'''
    SequenceDict = PassedSequenceDict #Separate from original dict to avoid pollution?
    assert isinstance(SequenceDict, dict)
    for key in ['title', 'sequence']:
        assert key in SequenceDict.keys()
    if 'accession' not in SequenceDict.keys():
        SequenceDict['accession'] = ''
    if 'comments' not in SequenceDict.keys():
        SequenceDict['comments'] = {}
    checkpoints = [0, len(SequenceDict['sequence'])]
    for point in SequenceDict['comments'].keys():
        assert isinstance(point, int)
        if point not in checkpoints:
            checkpoints.append(point)
    checkpoints.sort()
    Output = '>{0} {1}\r\n'.format(SequenceDict['accession'].strip(), SequenceDict['title'].strip())
    for point in checkpoints:
        if point != checkpoints[len(checkpoints)-1]: #i.e. if not last point
            nextpoint = checkpoints[checkpoints.index(point)+1]
        else: nextpoint = checkpoints[len(checkpoints)-1]
        SequenceFragment = SequenceDict['sequence'][point:nextpoint]
        FormattedSequence = ''
        for i in range(0, len(SequenceFragment))[::wraplines]:
            FormattedSequence += SequenceFragment[i:i+wraplines] + '\r\n'
        if point in SequenceDict['comments'].keys():
            Output += ';{0}\r\n{1}'.format(SequenceDict['comments'][point], FormattedSequence)
        else:
            Output += '{0}'.format(FormattedSequence)
    return Output.strip()
