# Sequence-Allignment-Dot-Plot
Uses Biopython and matplotlib to read in two Fasta sequences and create a sequence allignment dot plot with a window of ten. 

#Imports
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt

#Parsing through the two two files and assigning them to variables
seq1=list(SeqIO.parse(sys.argv[1],"fasta"))
seq2=list(SeqIO.parse(sys.argv[2],"fasta"))


#Getting the sequences and turning them into strings
#Also splices and stores the ids to be used as axis labes
for item in seq1:
    seq11=str(item.seq)
    xlab=item.id[10:]

for item in seq2:
    seq22=str(item.seq)
    ylab=item.id[10:]


#Creating two empty dictionaries
#Will store the subsets of the sequences and their indexes
dict_one={}
dict_two={}

#Setting the window size 
#window size: how many nucleotides you want to compare at one time
window=10

#looping through the sequence and appending 10 nucleotide long sections to dict_one
#sequence is the key and its indexes are the value 
for i in range(len(seq11)-window):
    section1=seq11[i:i+window]
    if section1 in dict_one:
        dict_one[section1].append(i)
    elif section1 not in dict_one:
        dict_one[section1]=[i]

#Same as above but for the second sequence
for i in range(len(seq22)-window):
    section2=seq22[i:i+window]
    if section2 in dict_two:
        dict_two[section2].append(i)
    elif section2 not in dict_two:
        dict_two[section2]=[i]


#contains the indexes of matches
xlst=[]
ylst=[]


#Loops through the keys and values of both dictionaires
#Creates a counter
#For every time the nucleotide from one key matches the nucleotide from the other key
#the count increases by one
#if the count is greater or equal to five the indexes are appended to the x or y lst    
for k,v in dict_one.items():
    for l in v:
        for k2,v2, in dict_two.items():
            for m in v2:
                count=0
                for i in range(0,len(k)):
                    if k[i]==k2[i]:
                        count=count+1            
                        if count>=5:
                            xlst.append(l)
                            ylst.append(m)
   



#sets up the scatter plot, dot size is set to two
plt.scatter(xlst, ylst, s=2)
#adds the axis labels
plt.xlabel(xlab)
plt.ylabel(ylab)
#adds the ticks (four ticks per axis)
plt.xticks(range(0,len(seq11),round(int(len(seq11))//int(4),-2)))
plt.yticks(range(0,len(seq22),round(int(len(seq22))//int(4),-2)))
#Needed to show the plot
plt.show()
