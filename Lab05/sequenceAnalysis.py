from FastAreader import FastAreader
from collections import OrderedDict
class NucParams:
    #Name: David Olmo Marchal

    """
    The NucParams class is the heart of preparing our genetic composition output, 
    It focuses on going through inputted genomes and making various different types of 
    compositions as like, amino acid composition, nucleotide composition, codon composition,
    and nucleotide count.
    
    """
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        '''
All the composition dictionaries are initialized in this constructor method of the NucParams class 
in which most are actually initialized to empty dictionaries 
to avoid wasting memory if you dont want the information that will go in it.        '''
        self.aacomposition = {value:0 for key,value in self.dnaCodonTable.items()}
        self.nucleotides = "" + inString
        self.codonComp = {}
        self.nucComp = {}
        # automatic call to addSequence since that is our way of updating the main running strand of code
        self.addSequence(self.nucleotides)
    def addSequence (self, inSeq):
        """
        Updates the current running code of genome assigned to the the NucParams object
        Checks for individual codons and seeing if its valid or not
        """
        self.nucleotides += inSeq
        if len(inSeq) >= 3:
            codon = ""
            #looping through bases
            for char in inSeq:
                codon += char
                #checking if valid codon
                if len(codon) == 3 and (codon in NucParams.dnaCodonTable.keys() or codon in NucParams.rnaCodonTable.keys()):
                    aa = NucParams.dnaCodonTable[codon]
                    #updating count of codon
                    self.aacomposition[aa] += 1 # inSeq.count(codon)
                    codon = ""
            #print(self.aacomposition)
        #updates nucleotide count
        for nuc in "ACGTNU":
            self.nucComp[nuc] = self.nucComp.get(nuc, 0) + inSeq.count(nuc)
        #print(self.nucComp)
    def aaComposition(self):
        '''
         returns already created amino acid composition
        '''
        return self.aacomposition 
    def nucComposition(self):
       '''
         returns already created and updated nuccomposition
        '''
       #print(self.nucComp)
       return self.nucComp
    def codonComposition(self):
        """specifically handles the count of specific codons inside the dna sequence"""
        #translates it to RNA coding
        RNAnucleotides = self.nucleotides.replace("T","U")
        if len(RNAnucleotides) >= 3:
            codon = ""
            #Loops through rna codons and counts
            for char in RNAnucleotides:
                codon += char
                if len(codon) == 3 and (codon in NucParams.rnaCodonTable.keys()):
                    self.codonComp[codon] = self.codonComp.get(codon,0) + 1 
                    codon = ""
        #print(self.codonComp)
        return self.codonComp
    def nucCount(self):
        ''' simple sum of the nucleotides'''
        return sum(self.nucComp.values())
    
class ORF:

    def __init__(self, startPos, stopPos, frame, length,sequence):
        ''' The point of this class is to contain different attributes of ORFs, like start/stop codon location, length of said ORF,
          just so it will be easy to sort later as we will have objects that already have everything we need'''
        self.startPos =startPos
        self.stopPos = stopPos
        self.frame = frame
        self.strand = frame[0]
        self.length = length
        self.sequence =sequence
    def __eq__(self, other):
        if not isinstance(other, ORF):
            return False
        return (self.startPos, self.stopPos, self.frame, self.length, self.sequence) == \
               (other.startPos, other.stopPos, other.frame, other.length, other.sequence)

    def __hash__(self):
        return hash((self.startPos, self.stopPos, self.frame, self.length, self.sequence))


class OrfFinder:

    def __init__(self,sequence,starts = ['ATG','GTG','CTG','TTG'],stops=['TAA','TAG','TGA']):
        '''In this constructor we focus on preparing the different DNA sequences we will need to look at for our ORFs.
        This includes the normal forward strand, the reverse complement strand and the 3 different frame shifts of those two strands
        So a total of 6 different strands we are looking at for the ORFs'''
        self.sequence = sequence
        self.starts = starts
        self.stops = stops
        self.orfDict = {}
        #print("This is the normal sequence: %s"%sequence)
        #We need to find the complement, so we need to translate bases
        complementmap = str.maketrans('ACGT', 'TGCA')
        self.complement = sequence.translate(complementmap)
        #print("This is the complement sequence: %s"%self.complement)
        # turning it into a list so we can use reverse method
        complementlist = list(self.complement)
        complementlist.reverse()
        self.revcomplement = "".join(complementlist)
        #print("This is the reverse complement sequence: %s"%self.revcomplement)
        self.frames = {
          '+1':self.frameShift(self.sequence,1),
          '+2':self.frameShift(self.sequence,2),
          '+3':self.frameShift(self.sequence,3),
          '-1':self.frameShift(self.revcomplement,1),
          '-2':self.frameShift(self.revcomplement,2),
          '-3':self.frameShift(self.revcomplement,3)
        }
        #print("This is the frame +1 : %s"%self.frames.get("+1"))
        #print("This is the frame +2 : %s"%self.frames.get("+2"))
        #print("This is the frame +3 : %s"%self.frames.get("+3"))
        #print("This is the frame -1 : %s"%self.frames.get("-1"))
        #print("This is the frame -2 : %s"%self.frames.get("-2"))
        #print("This is the frame -3 : %s"%self.frames.get("-3"))

    def frameShift(self,sequence,frame):
        '''This helper function changes the frame of the given sequence
        depending on what has been given by as an input
        Ex. Frame +1 ATG CCT TAG GCT GA  
            Frame +2 TGC CTT AGG CTG A
            Frame +3 GCC TTA GGC TGA'''
        #shifts the frame
        frame = frame - 1
        return sequence[frame:]
    
    def isStart(self, codon):
        '''checks if its a start codon option'''
        return  (codon in self.starts)
    def isStop(self,codon):
        '''checks if its in the stop codon options'''
        return (codon in self.stops)
    def sortOrfs(self):
        orfs = list(self.orfDict.items())
        orfs.sort(key=lambda x: (-x[1].length, x[1].startPos,x[1].frame))  # Sort by length (desc), then startPos (asc)
        self.orfDict = OrderedDict(orfs)
    def remove_duplicate_orfs(self):
        unique_orfs = set(self.orfDict.values()) 
        new_orf_dict = {i: orf for i, orf in enumerate(unique_orfs)} 
        self.orfDict = new_orf_dict
    def orfFinder(self, frame,mG = 0):
        '''Start codons are ATG|GTG|CTG|TTG 
            Stop codons are TAG|TAA|TGA '''
        # there will be a dictionaries filled with ORF objects which 
        # we will be easily be able to access, sorted by length
        #lists that handle the positions of start and stop codons
        # gets the specific frame sequence from dictionaries
        
        frameSeq = self.frames.get(frame)
        beginning = True
        foundStartorStop = False
        for i in range(0,len(frameSeq)-2,3):
            #loops through each odon and then checks to see if its start or stop
            codon = frameSeq[i:i+3]
            #boundary condition
            if self.isStop(codon) and beginning == True:
                # This is the test case for when we find a stop codon in the beginning but no start codon
                # This means that we have an ORF from the start of the sequence until that stop codon (assume start is upstream)
                stopPos = i+3
                length = stopPos
                if(length>=mG):
                    if frame[0] == "+":
                        boundarySeq = self.frames.get("+1")
                    else:
                        boundarySeq = self.frames.get("-1")
                    seq = boundarySeq[0:stopPos+(int(frame[1])-1)]
                    length = len(seq)
                    nameOrf= "ORF "+ str(len(self.orfDict)+1)
                    orf = ORF(1,length,frame,length,seq)
                    self.orfDict[nameOrf] = orf
                    foundStartorStop = True
                beginning = False
            
            if self.isStart(codon):
                #beginning = False
                #This is for finding ORFs in the middle of the sequence, just like ur average ORF finder, so no boundary conditions
                startPos = i
                for j in range(i, len(frameSeq)-2,3):
                    stop_codon = frameSeq[j:j+3]
                    if self.isStop(stop_codon):
                        #if a stop codon append it to the list
                        stopPos = j+3
                        seq = frameSeq[startPos:stopPos]
                        length = len(seq)
                        #creates a ORF object and adds it to our dictionaries of ORFs read
                        if(length>=mG):
                            nameOrf= "ORF "+ str(len(self.orfDict)+1)
                            orf = ORF(startPos+1,stopPos,frame,length,seq)      #fix this later might be hard coded
                            self.orfDict[nameOrf] = orf
                            foundStartorStop = True
                            break
        if foundStartorStop != True:
            #This is the boundary condition if u find zero start or zero stop codons
            if frame[0] == "+":
                boundarySeq = self.frames.get("+1")
            else:
                boundarySeq = self.frames.get("-1")
            length = len(boundarySeq)
            if(length>=mG):
                nameOrf= "ORF "+ str(len(self.orfDict)+1)
                orf = ORF(1,length,frame,length,boundarySeq)
                self.orfDict[nameOrf] = orf
                foundStartorStop = True
                beginning = False
            
        return self.orfDict

def main():

    myReader = FastAreader('lab5test.fa') 
    for head, seq in myReader.readFasta() :
        print(head)
        orf = OrfFinder(seq)
        orf.orfFinder("+1")
        orf.orfFinder("+2")
        orf.orfFinder("+3")
        orf.orfFinder("-1")
        orf.orfFinder("-2")
        orf.orfFinder("-3")
        orf.remove_duplicate_orfs()
        orf.sortOrfs()

        
        for orf_name, orf in orf.orfDict.items():
            print(f"{orf.frame:<4} {orf.startPos:>3}..{orf.stopPos:>3} {orf.length:>4} {orf.sequence}")
        
            #for key,value in orftest.orfDict.items():
                #print("%s: Strand: %s Frame: %s Start %s Stop %d Length(nt|aa) %d | %d"%(key, value.strand, value.frame, (value.startPos ), (value.stopPos ), value.length, (value.length/3)))
    
if __name__ == "__main__":
    main()