import sys

# This script allows to remove all N character in a sequence file and replace it with a gap (-)

in1 = sys.argv[1] # this is the reading file
out = sys.argv[2] # this is the output file

with open(in1, "r") as infile, open(out, "w") as outfile:
    lines = infile.readlines()
    for line in lines :
        if line.startswith(">"):
            outfile.write(line)
        else:
            for elt in line:
                if elt == "n" or elt == "N":
                    outfile.write("-")
                else:
                    outfile.write(elt)
