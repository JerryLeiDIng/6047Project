
#Set up the list of SNPs we need
f = open("GWAS.txt", 'r')

gwasSnps = []

line = f.readline()
while (line != ""):
    line = f.readline()
    tokens = line.split(" ")
    gwasSnps.append((tokens[0],tokens[2]))

f.close()

ld = f.open("LD.txt", 'r')

ldcoeffs = {}
allSnps = []

line = ld.readline()
currentSnp = None
while(line != ""):
    tokens = line.split(" ")
    if not tokens[3] == currentSnp:
        currentSnp = tokens[3]
        allSnps.add(currentSnp)
    ldcoeffs[(currentSnp,tokens[4])] = (tokens[6],tokens[1]) #TODO
    line = ld.readline()
    
ld.close()

snpRanges = []
#Builds snpRanges: list of tuples (snp, minPos, maxPos of linked alleles)
for snpPair in gwasSnps:
    snp = snpPair[0]
    maxPos = snpPair[1]
    minPos = snpPair[1]
    lds = ldcoeffs[snp]
    for otherSnp in allSnps:
        if otherSnp == snp:
            continue
        else if (otherSnp, snp) in ldcoeffs.keys():
            current = ldcoeffs[(otherSnp,snp)]
        else:
            current = ldcoeffs[(snp,otherSnp)]
            
        if current[0] >= 0.5:
            if current[1] > maxPos:
                maxPos = current[1]
            if current[1] < minPos:
                minPos = current[1]

    snpRanges.append((snp,minPos,maxPos))

#Write the SNP ranges to a file. Can get rid of this later
snpRangeFile = open("snp_ranges",'w')
for snp in snpRanges:
    snpRangeFile.write("{0[0]} {0[1]} {0[2]}".format(snp))
snpRangeFile.close()

#TODO associate genes with SNPs and write to file

