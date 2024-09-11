#Author: Moritz berger
#Date: 06.09.2023

#imports: 
import csv
import sys

yseqid = sys.argv[1]


#Script to create YSEQID_mtDNA.fasta

#read mtdna fasta file and split into header and dna.
#also remove all newline [\n] from the sequence

fasta_file = open("/genomes/0/refseq/mt/rCRS.fa", "r")
#fasta_file = open("test-ref.fasta", "r")

header = fasta_file.readline()

# Change header (We don't want to keep the rCRS header)
header = ">chrM " + str(yseqid) + "\n"

dna = fasta_file.read().replace("\n","")


#read YSEQID_MTDNA_SNPS.tsv
snps_f=open(str(yseqid)+"_MTDNA_SNPS.tsv")

#snps_f=open("test.tsv")
tsv_reader = csv.reader(snps_f, delimiter='\t')


#skip header line in tsv
next(tsv_reader, None)

#for each pos & base, add a Tupel to the List
mutation = []
for row in tsv_reader:
	pos, base = row
	base = (base.upper().strip())[:1]
	if (pos.find('.')!=-1) :
		pos = pos.split('.')[0]
		base = base +'+'

	mutation.append((int(pos)-1, base))

for pos,base in mutation:
	print("position: "+str(pos) + " Base: " + str(base) )

#reverse iterate the List of Tupels 
for pos, base in reversed(mutation):
	#print(pos,"\t",base)

	#If deletion[D] then delete the pos.
	if (base in ["D", "DEL", "d", "del"]):
		dna = dna[:pos]+dna[pos +1:]
	#If insertion [X+], insert base at after pos
	elif (base.find('+') == 1):
		dna = dna[:(pos+1)] + base.split('+')[0]+ dna[pos+1:]
	#else replace at pos with 
	else:
		dna = dna[:pos] + base + dna[pos+1:]


#split string in lines of 70 characters
char ="\n"
group = 70
dna = char.join(dna[i:i+group] for i in range(0, len(dna), group))

#save fasta "yseqid_mtDNA.fasta"
out_file = open(str(yseqid)+"_mtDNA.fasta", "w")
out_file.write(header)
out_file.write(dna)
