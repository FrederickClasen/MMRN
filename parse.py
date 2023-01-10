import sys

fin = open(sys.argv[1],'r')
fout = open(sys.argv[2],'w')

for line in fin:
    if '<species metaid="MMRNM' in line:
        parts = line.split()
        comp = str(parts[1][-2]).strip()
        line = line.replace('compartment="c"','compartment='+'"'+comp+'"')
        fout.write(line)
    else:
        fout.write(line)
fin.close()
fout.close()

