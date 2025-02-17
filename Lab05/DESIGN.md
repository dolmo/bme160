## Brainstorming
    - I want to use regular expressions, just cause I feel like they will be way more efficient the given the size of the genome we are given.
    - But if regular expressions dont work I will move over to indexing. (Post-planning edit. Ended up taking this route because of a lack of time and difficulty of debugging regex, will include the brainstorming for itthough)
# Methods I think I should use
    - re.search()
    - .start()
    -re.split()
# Some regular expressions that might  be useful
    - r"([ACTG])+ (TAG|TAA|TGA)"
        - This is going to match a DNA Sequence terminated by a stop codon. ([ACTG])+ will get all but the last stop codon.
        - When we put this into https://regex101.com/ and see how it reacts to some genome we see that it splits the genome into groups whnever it reaches the stop codon multiple times
        - We might want to use this version of the regex:
        - r"([ACTG])+? (TAG|TAA|TGA)"
        - Since this splits once it finds ONE stop codon.
    - r"([ACTG])+? (ATG|GTG|CTG|TTG)"
    - Stripping white space
        -  re.split(r"\s+",line.strip())
 # Sources
    - https://open.oregonstate.education/computationalbiology/chapter/bioinformatics-knick-knacks-and-regular-expressions/
    - https://regex101.com/

# Pseaudocode
1. Get the DNA String ready for use. *
2. Get the reverse complement, this involves going from A->T, G->C and C->G, then you need to reverse the order *
3. Now that you have two strands you need to look at the different frames
4. There are 3 Frames for each strand (forward and reverse complement), you need to just move up by +1 nucleotide for a diff frame shift. Ex ATGCCTTAGGCTGA, Frame +1 → ATG CCT TAG GCT GA, Frame +2 → TGC CTT AGG CTG A, Frame +3 → GCC TTA GGC TGA., you do this with reverse aswell.
5. Now we find the start and end codons, and begin the search for specific ORFs.
6. Im going to make an ORF class, which will contain different attributes of ORFs, like start/stop codon location, length of said ORF, just so it will be easy to sort later as we will have objects that already have everything we need.

# Things to think about
1. We can take into account a minimum gene size, for example 100 nucleotides (counts start and stop codon so 94 coding Nuc and 6 start/stop nuc)
2. We can pick out the largest ORF, just simply taking the size of the ORF. So we make a method that figures out the max. Like the C reduce method.
3. Have variable start codons- this is easy
4. Have variable stop codons - this is easy
5. Using indexing we can easily calculate the new frame shifts, (position % 3) + 1.