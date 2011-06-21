from optparse import OptionParser, OptionGroup
from Bio import SeqIO
import DDB
import sys
usage = "usage: %prog fasta_file outputfile\nAfterwards run SSRcalc:\n" 
usage += "perl external/SSRCalc4.pl --alg 3.0 --source_file outfile  --output tsv --B 1 --A 0  >> ssrcalc.out"
parser = OptionParser(usage=usage)
options, args = parser.parse_args(sys.argv[1:])

fasta_file = sys.argv[1]
outfile = sys.argv[2]
records = list(SeqIO.parse(open(fasta_file,"r"), "fasta"))

done_already = {}
f = open(outfile, 'w')
for r in records:
    peptides = DDB.Protein.trypsin(r.seq.tostring() )
    for p in peptides: 
        if not done_already.has_key(p): f.write('%s\n' % p); done_already[p] = 0

f.close()

