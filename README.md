# PySplicer
by Cathal Garvey, August 2012.
Documentation and JSON under Creative Commons, Attribution, Sharealike License.
Python code under GNU Public License.

## Que?
PySplicer is a codon optimisation script that operates on a weighted random basis to try and provide an RNA sequence for a given amino acid based on a desired codon frequency table.

When used with a codon frequency table for constitutively expressed genes in a target species, this is a fairly good way of codon-optimising a sequence for expression in a target species. When used with a table that has been empirically proven to provide high expression, this can be used to design genes for protein overexpression and production.

I wrote this because there are not many tools that support this approach; most codon optimisation tools available today are based on a somewhat disproven "best pick" approach that preferentially chooses the "most adapted" codon for each amino acid. Experimental evidence shows that this approach not only fails to achieve exceptional results, but can sometimes lead to *worse* protein expression, perhaps due to exhaustion of a limited pool of charged tRNAs for a given codon faster than cognate tRNA synthetases can compensate.

## Usage
PySplicer is a CLI application as written, although the code should be adaptable to web services using something like CherryPy without too much extra difficulty.
CLI usage can be informed by typing "./pysplicer.py -h" in terminal, when in the same directory as pysplicer.

Suggested usage is something like this:
me@mycomputer:~/Applications/PySplicer$ ./pysplicer.py -x Enzymes_Example GFP_WT

If you want to see *loads* of debugging information, try adding the -v argument.

## Process
PySplicer first generates a set of candidate RNA codon lists from the given amino acid sequence. By default, it creates 20 sets, but you can change this with the -n argument in the Terminal.

Each of these lists is created by selecting by weighted-random a codon for each amino acid. You can choose a minimum frequency, whereupon the optimiser will recalculate the codon table as if any codons less common than the given frequency simply don't exist. This may cause crashes, see "unresolved issues" below.

The optimiser then maps all the undesirable enzyme sites in each candidate and attempts to "splice" a hybrid sequence lacking all of these sites. When it is finished (hopefully with no remaining sites), it then performs a sanity check to ensure the codons still encode the specified amino sequence, and passes the candidate onto a dedicated site excision method.

In this method, the candidate is again searched for target enzyme sites; these might remain from splicing, or may have been created accidentally at splice junctions. This method attempts to resolve enzyme sites by generating permuations for each detected site, and randomly selecting a permutation that doesn't introduce new sites into the codon list.

At present, selection of replacement codons is entirely random; a planned feature is a frequency evaluation and comparison to the desired codon frequency table, with preference being shown to replacement codons that bring the candidate codon list closer to the desired frequency spread.

## Advantages
* This script is Free/Libre Open Source Software, and appears to be the first FLOSS script to address this method of codon optimisation.
* This approach allows you to rapidly exclude degenerate sites while optimising an RNA sequence for a given amino acid sequence.
* This script can parse degenerate IUPAC sequences for sites to exclude (which are given as DNA sequences, rather than RNA).
* This script comes with an enzyme profile picker which (while very dumb) makes it easy to add sites for most known restriction enzymes from an extensive master list.
* With the right frequency tables, it is thought that this approach can offer the best translation efficiencies currently achievable by publicly disclosed means.
* This object/script is portable to other applications, for example mobile apps though Android Scripting Layer, webapps through CherryPy, or chat bots through python programs such as Phenny.

## Disadvantages
* This script only supports FASTA files at present.
* This script currently has no regard for the original translation speed of the wild gene. It is now known that some proteins require "poorly adapted" codons at certain locations in order to slow down translation, allowing proper folding of the protein structure before further amino acid addition. With this script, no such adaptations will be respected, and some rare proteins may therefore suffer.
* This script has no frequency analysis routines at present, and so it is left to the user to ensure that statistical fluke does not lead to a terribly "optimised" sequence. Frequency analysis is a desired feature in future revisions of the code. At the very least, a verbose-mode printout at the end listing desired codon frequencies versus delivered frequencies would be helpful.
* This script requires codon tables to be provided in a custom-build JSON format. While very simple, this necessitates converting a desired codon table in one of the many hard-to-parse formats in widespread use into the JSON format used by this script. Also, this script uses codon frequencies relative to synonymous codons (that is, if an amino acid is encoded by two codons, and they are used in an 80:20 ratio, the frequencies will be 0.8 and 0.2, respectively), whereas it is customary to use frequencies per 1000 codons in most frequency tables. This will also necessitate conversion. A helper script to automate this tedious process is another desired feature.
* The default codon frequency table is derived from evidence based work, but only in part; the study in question didn't provide data on all codons, and so the "gaps" were filled in with the genome-wide codon frequencies for *E.coli K12*. This may be a bad decision, and perhaps someone knows better how to fill gaps in data such as this for optimal expression.

## Unresolved Issues and Desired Features
* Setting a minimum codon frequency above a certain threshold, which likely varies from table to table, will cause the script to fail.
* A refactoring of the "map_excluded_sites" method to use regular expressions rather than straight lists of comparison substrings would avoid issues where the program might crash if given *highly* degenerate sequences like AGGANNNNNNNNNNWWWWWG, which would compute to no less than 33,554,432 substrings. The above could be represented by a regex as "AGGA[ACGT]{10}[AT]{5}G" and would require a lower memory footprint and likely less time to parse.

## Credit Where It's Due
* The master list of enzymes provided with the script for ease of exclusion profile creation is derived from NEB's REBASE database of Restriction Enzymes. It's provided in JSON format and may prove useful to others as such, but it wouldn't have been possible to compile without REBASE itself.
* The empirical work that informed the creation of the default frequency table for *E.coli* is courtesy of Welch et al, 2009, who created libraries of artificial genes and used them to verify that expression was best predicted by using codons whose tRNAs remained charged under starvation conditions, rather than the most common codons in "highly expressed" or "constitutive" genes. I look forward to seeing similar work done in other common target species.
