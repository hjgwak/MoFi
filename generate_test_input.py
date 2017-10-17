import sys, random

def getArg(argv, key, default) :
    if key in argv :
        return argv[argv.index(key) + 1]
    return default

def makeMotifProfile(k) :
#    template = [{'A' : 0.8,  'C' : 0.1,  'G' : 0.05, 'T' : 0.05},\
#                {'A' : 0.1,  'C' : 0.8,  'G' : 0.05, 'T' : 0.05},\
#                {'A' : 0.05, 'C' : 0.05, 'G' : 0.8,  'T' : 0.1},\
#                {'A' : 0.05, 'C' : 0.05, 'G' : 0.1,  'T' : 0.8},\
#                {'A' : 0.5,  'C' : 0.2,  'G' : 0.2,  'T' : 0.1},\
#                {'A' : 0.2,  'C' : 0.1,  'G' : 0.4,  'T' : 0.3},\
#                {'A' : 0.15, 'C' : 0.25, 'G' : 0.25, 'T' : 0.35},\
#                {'A' : 0.1,  'C' : 0.6,  'G' : 0.2,  'T' : 0.1},\
#                {'A' : 0.1,  'C' : 0.1,  'G' : 0.2,  'T' : 0.6}]
    template = [{'A' : 1.0, 'C' : 0.0, 'G' : 0.0, 'T' : 0.0},\
                {'A' : 0.0, 'C' : 1.0, 'G' : 0.0, 'T' : 0.0},\
                {'A' : 0.0, 'C' : 0.0, 'G' : 1.0, 'T' : 0.0},\
                {'A' : 0.0, 'C' : 0.0, 'G' : 0.0, 'T' : 1.0}]
    profile = []
    for i in range(k) :
        profile.append(template[random.randint(0, len(template) - 1)])

    return profile

def makeMotif(profile) :
    motif = ""
    for i in range(len(profile)) :
        die = random.randint(0, 10000) / 10000.0
        if die <= profile[i]['A'] :
            motif += 'A'
        elif die <= profile[i]['A'] + profile[i]['C'] :
            motif += 'C'
        elif die <= profile[i]['A'] + profile[i]['C'] + profile[i]['G'] :
            motif += 'G'
        else :
            motif += 'T'
    return motif

# main

if len(sys.argv) < 3 :
    print "python make_test_input.py [-k #] [-l #] [-n #] [-o str]"
    print "       -k\t:\tmotif size, default 10"
    print "       -l\t:\tsequence length, default 100"
    print "       -n\t:\tnumber of sequences, default 50"
    print "       -o\t:\toutput file name, required"
    exit(-1)

k = int(getArg(sys.argv, "-k", "10"))
l = int(getArg(sys.argv, "-l", "100"))
n = int(getArg(sys.argv, '-n', "50"))
o = getArg(sys.argv, "-o", "")

if o == "" :
    print "output file name is requried"
    exit(-1)

if l < k :
    print "Motif size should be shorter than sequence length"
    exut(-1)

print k
print l
print n
print o
output = open(o, 'w')
motif_ = open(o + ".motif", 'w')

profile = makeMotifProfile(k)

for i in range(n) :
    motif = makeMotif(profile)
    motif_pos = random.randint(0, l - k)
    
    output.write(">test_" + str(i) + "\n")
    motif_.write(">test_" + str(i) + "\n")
    for j in range(motif_pos) :
        output.write(random.choice(['A', 'C', 'G', 'T']))
    output.write(motif)
    motif_.write(motif + "\n")
    for j in range(l - motif_pos - k) :
        output.write(random.choice(['A', 'C', 'G', 'T']))
    output.write('\n')

output.close()

output = open(o + ".profile", 'w')
for ch in ['A', 'C', 'G', 'T'] :
    out = []
    for i in range(k) :
        out.append(str(profile[i][ch]))
    output.write(",".join(out) + "\n")
output.close()
