
# This script allows to remove all IUPAC character in a sequence file and replace it with a gap (-)

import sys

in1 = sys.argv[1] # this is the reading file
out = sys.argv[2] # this is the output file
nc = ["A", "T", "C", "G"]

with open(in1, "r") as infile, open(out, "w") as outfile:
    lines = infile.readlines()
    for line in lines :
        line.upper()
        if line.startswith(">"):
            outfile.write(line)
        else:
            for elt in line:
                if elt not in nc:
                    outfile.write("-")
                else:
                    outfile.write(elt)

