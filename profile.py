import sys

f = open(sys.argv[-1], 'r')

profile = []
for i in range(int(sys.argv[-2])) :
    profile.append({'A' : 0.0, 'C' : 0.0, 'G' : 0.0, 'T' : 0.0})

total = 0.0
for line in f.readlines() :
    line = line.rstrip('\r\n')
    if line[0] == '>' :
        continue
    for i in range(len(line)) :
        profile[i][line[i]] += 1.0

    total += 1.0

f.close()

for prob in profile :
    prob['A'] /= total
    prob['C'] /= total
    prob['G'] /= total
    prob['T'] /= total

for ch in ['A', 'C', 'G', 'T'] :
    out = []
    for prob in profile :
        out.append(str(prob[ch]))
    print ",".join(out)
