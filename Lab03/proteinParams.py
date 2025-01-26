# Name: David Olmo Marchal

class ProteinParam :
    """
    In this class we create a protein object that has multiple methods to analyze its composition and the different chemical properties it has
        -We count the amount of valid amino acids
        -Make a dictionary of the different amino acids with amount values attached
        -We calcualte the net charge of the protein at a certain pH
        -We calculate what theoritical pH would make our protein have a charge closest to 0.0
        -Calculate other misc chemical properties like molar extinction, mass extinction and molecular weight
    
    """
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        
        """Initializing the different instance variables we will be using throughout the class"""
        self.protein = protein
        self.aa2mw = ProteinParam.aa2mw
        self.pkaPos = ProteinParam.aa2chargePos
        self.pkaNeg = ProteinParam.aa2chargeNeg

        #Using dictionary comprehension we make a count of the number of specific amino acids inside our protein  and put  it in a dictionary
        self.proteincomposition = {aa: protein.count(aa) for aa in self.aa2mw.keys()}
    #here we count the number of valid amino acids in our protein and remove invalid amino acids
    def aaCount (self):
        """Gets the count of total valid amnino acids in the protein"""
        #looping through the protein string
        for item in self.protein:
            #checking if its invalid or not, if its not in the keys for our amino acid molecular weight then its not a valid amino acid
            if item not in self.aa2mw.keys():
                self.protein.replace(item, "")
        #returns the shortened length of the list of amino acids in our protein
        return len(self.protein)
    

    def aaComposition (self) :
        """simply returns the composition we made earlier in the __init__ method"""
        return self.proteincomposition

    def _charge_ (self,pH):
        #Defining important variables for calculation
        positive = 0.0
        negative = 0.0
        netCharge = 0.0
        #Looping throught the proteincomposition dictionary yet again
        for key,value in self.proteincomposition.items():
            #checking for eligible amino acids for pka in the positive range
            if value > 0 and key in self.pkaPos.keys():
               #calculating the positive charge
               positive += float(value) * 10**self.pkaPos.get(key) / (10**(self.pkaPos.get(key))+ 10**pH)
            #checking for eligible amino acids for pka in the negative range
            elif value > 0 and key in self.pkaNeg.keys():
                #calculating the negative charge 
                negative += float(value) * 10**pH/ (10**(self.pkaNeg.get(key))+ 10**pH)
        #Nterminus and Cterminus taken into account
        positive += 10**ProteinParam.aaNterm /(10**(ProteinParam.aaNterm)+ 10**pH)
        negative += 10**pH /(10**(ProteinParam.aaCterm)+ 10**pH)

        netCharge = positive-negative
        return netCharge

    
    def pI (self):
        """In this method we calculate theoritical pH needed to get a net charge around 0.0
            for this method I had tried doing binary search but I had a lot problems with it so I switched to a linear approach """
        pH=0.0
        nearest0 = float('inf')
        bestpH = 0.0
        #Looping through the 0-14 pH range
        while pH <=14.0:
            #calling the charge method given the current pH index
            charge = self._charge_(pH)
            #checking if the current charge is less than the previous lowest charge
            if abs(charge) < nearest0:
                # we use abs because netcharge can mean a negative charge and could mess up the comparisons
                nearest0 = abs(charge)
                #updates the bestpH variable everytime we get a new closer pH
                bestpH = pH
            #up our pH by 0.1 everytime
            pH += 0.01
        return bestpH
            
    def molarExtinction (self):
        """in this method  we calculate the molar extinction coefficient by looping through our protein composition and following the formula given"""
        me = 0.0
        #looping through our protein composition dictionary again
        for key,value in self.proteincomposition.items():
            #check for amount of the 3 specific amino acids
            if value > 0 and key in self.aa2abs280.keys():
                #calculation
                me += value *(self.aa2abs280.get(key))
        return me

    def massExtinction (self):
        """not written by me but it calculates massExtinction by calling molecularWeight method."""
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        """In this method we calculate the molecular weight of the counted amino acids """

        mw = 0.0
        #Looping through the dictionary based on keys and values 
        for key,value in self.proteincomposition.items():
            #Every amino acid that has at least one of itself in the protein
            if value > 0:
                #those amino acids with at least one are substracted by the molecular weight of water and then multiplied by the numebr of amino acids they have
                mw += value *(self.aa2mw.get(key)-ProteinParam.mwH2O) 
            #Then we return the value fo that summed molecular weight plus one included water molecular weight
        return ProteinParam.mwH2O + mw 
# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    #remove later
    inString = input('protein sequence?')
    myParamMaker = ProteinParam(inString)
    myAAnumber = myParamMaker.aaCount()
    
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()