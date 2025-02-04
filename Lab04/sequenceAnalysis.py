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