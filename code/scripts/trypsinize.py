"""
Input:

    Fasta file format:

>2977010
MAPVVISESEEDEDRVAITRRTKRQVHFDGEGDDRVDQQQQQHSSSHRDRDKHVQRKKKKRLSNRNLQGSNGGYAWEDEIK
RSWDLVKVDDEGDMASLVASIVEARKKRTAKKNITPYQRGIIRSLILTLDCSEAMLEKDLRPNRHAMIIQYAIDFVHEFFD
QNPISQMGIIIMRNGLAQLVSQVSGNPQDHIDALKSIRKQEPKGNPSLQNALEMARGLLLPVPAHCTREVLIVFGSLSTTD
>2975962
MLRNGVQRLYS...

Output:
MAPVVISESEEDEDR
MAPVVISESEEDEDRVAITR
VAITR
VAITRR
RQVHFDGEGDDR
SWDLVK
...

"""
from optparse import OptionParser, OptionGroup
import sys; sys.path.extend(['..', '.'])
from Bio import SeqIO

usage = 'A script to read a fasta file and output trypsinized peptides, one per line\n'
usage += "usage: %prog fasta_file outputfile missed_cleavages min_len\nAfterwards run SSRcalc:\n" 
usage += "perl SSRCalc3.pl --alg 3.0 --source_file peptide_file  --output tsv --B 1 --A 0  > ssrcalc.out"
parser = OptionParser(usage=usage)
options, args = parser.parse_args(sys.argv[1:])

fasta_file = sys.argv[1]
outfile = sys.argv[2]
missed = 0
min_len = 0
if(len(sys.argv)>3):
    missed = int(sys.argv[3])

if(len(sys.argv)>4):
    min_len = int(sys.argv[4])

records = list(SeqIO.parse(open(fasta_file,"r"), "fasta"))

def trypsinize(sequence, missed=0):
    import re
    protein = re.sub(r'(?<=[RK])(?=[^P])',' ', sequence)
    protein = protein.split()
    for i,peptide in enumerate(protein):
      yield peptide 
      # do missed cleavages
      current = peptide
      k = 1
      while(i+k<len(protein) and k<=missed):
          current = current + protein[i+k]
          yield current
          k+=1

done_already = {}
f = open(outfile, 'w')
for r in records:
    peptides = trypsinize(r.seq.tostring(), missed)
    for p in peptides: 
        if not done_already.has_key(p) and len(p) >= min_len:
            f.write('%s\n' % p); done_already[p] = 0

f.close()

"""
mysql -s -e "select protein.id,sequence from ddbMeta.sequence inner join \
        ddb.protein on sequence.id = protein.sequence_key where experiment_key \
        = 3131 limit 2" | perl -ane 'printf ">%d\n%s\n", @F; '

"""
