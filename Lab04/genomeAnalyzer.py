#Name: David Olmo Marchal
from  sequenceAnalysis import NucParams
from FastAreader import FastAreader


def main (fileName=None):
    '''
    This is where we handle all the output that we get from fasta files like the testgenome.fa file
    We use the nucparams class and format it into specific values we want to see
    like amount of Mb, GC compositon or codon composition.
    '''
    #objects created
    myReader = FastAreader(fileName) 
    myNuc = NucParams()
    #looping through file and calling add sequence
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)
    #some placeholder variables
    codons = myNuc.codonComposition()
    nucleotideCount = myNuc.nucCount()
    #printing sequence length with right sig figs
    print("sequence length = %.2f"%(nucleotideCount/1000000) + " Mb")
    nuccomp = myNuc.nucComposition()
    #gc composition
    gccontent = ((nuccomp.get("G") + nuccomp.get("C")) / nucleotideCount) * 100
    print("\nGC content = %.1f"%(gccontent) + "%\n")
    # sort codons in alpha order, by Amino Acid
    
    # calculate relative codon usage for each codon and print
    #orders
    for codon,count in sorted(codons.items(), key=lambda x: (NucParams.rnaCodonTable.get(x[0]),x[0])):
        aa = NucParams.rnaCodonTable.get(codon)
        total = sum(codons[c] for c in codons if NucParams.rnaCodonTable.get(c) == aa)
        codonval = (count / total)
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, NucParams.rnaCodonTable.get(codon), codonval*100, count))
if __name__ == "__main__":
    main('testGenome.fa') # make sure to change this in order to use stdin
    