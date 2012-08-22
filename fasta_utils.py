#!/usr/bin/env python3
'''Some functions for working with FASTA sequences as native types.'''

def read_fasta_file(file):
    'Imports the contents of a file as a string, passes to read_fasta() and returns output.'
    with open(file) as fasta_file:
        fasta_dict = read_fasta(fasta_file.read())
    return fasta_dict

def read_fasta(fasta_string):
    '''Finds a sequence and returns a dictionary containing title, comments, accession (if found) and sequence.

    Supports extraction of titles, comments and sequence. (Comments are lines beginning with ';'.
    This is an oft-ignored FASTA spec feature and is almost never used.)
    Dictionary returned contains following keys, whether empty (evals false) or with contents:
     'title': (a string containing the title, including the ">" symbol)
     'comments': (a dict of format {string index: comment})
     'sequence': (a string containing the sequence)
     'accession': (a string containing the accession information if present (detected by searching for '|' in title)
    '''
    sequence = {'title':'','comments':{},'sequence':'','accession':''}
    fasta_lines = []
    for line in fasta_string.split("\n"):
        line = line.strip()
        if line.startswith(">"): # Is sequence descriptor
            sequence['title'] = line[1:].strip()
        elif line.startswith(";"): # Is a comment, rarely used
            sequence['comments'][len(sequence['sequence'])] = line[1:].strip()
            # Uses sequence length SO FAR minus 1, so acts like a string index for comment markup location.
            # sequence['comments'].append(Comment)
        else:
            sequence['sequence'] += line.upper().strip()
    first_title_word = sequence['title'].split(' ')[0].strip()
    if '|' in first_title_word: # Yes, this is crude, I know. #lazy
        sequence['accession'] = first_title_word
        sequence['title'] = sequence['title'][len(sequence['accession']):].strip()
    return sequence

def make_fasta(passed_sequence_dict, wraplines = 50):
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
    sequence_dict = Passedsequence_dict #Separate from original dict to avoid pollution?
    assert isinstance(sequence_dict, dict)
    for key in ['title', 'sequence']:
        assert key in sequence_dict.keys()
    if 'accession' not in sequence_dict.keys():
        sequence_dict['accession'] = ''
    if 'comments' not in sequence_dict.keys():
        sequence_dict['comments'] = {}
    checkpoints = [0, len(sequence_dict['sequence'])]
    for point in sequence_dict['comments'].keys():
        assert isinstance(point, int)
        if point not in checkpoints:
            checkpoints.append(point)
    checkpoints.sort()
    output = '>{0} {1}\r\n'.format(sequence_dict['accession'].strip(), sequence_dict['title'].strip())
    for point in checkpoints:
        if point != checkpoints[len(checkpoints)-1]: #i.e. if not last point
            nextpoint = checkpoints[checkpoints.index(point)+1]
        else: nextpoint = checkpoints[len(checkpoints)-1]
        sequence_fragment = sequence_dict['sequence'][point:nextpoint]
        formatted_sequence = ''
        for i in range(0, len(sequence_fragment))[::wraplines]:
            formatted_sequence += sequence_fragment[i:i+wraplines] + '\r\n'
        if point in sequence_dict['comments'].keys():
            output += ';{0}\r\n{1}'.format(sequence_dict['comments'][point], formatted_sequence)
        else:
            output += '{0}'.format(formatted_sequence)
    return output.strip()
