#put the DNA sequence from the 5- direction and without 5- or 3-:
DNA= input('Please put your sequence.').upper().strip() #DNA sequence input
x =(f"The Main Strand is: \n 5- {DNA} -3\n") #DNA strand output 
def complementry(DNA): 
    "give the complementry strand of the DNA"
    trans = DNA.maketrans("ATGC", "TACG") # maketrans make a map table and use it in translate() to to replace specified characters.
    global strand #use global variable to use it out the funcation 
    strand = f'{DNA.translate(trans)}'.strip() #the complementry DNA 
    return f'The Complementry Strand is: \n 3- {strand} -5\n'
def transcription(DNA):
    "convert DNA to to RNA"
    trans = DNA.maketrans("GCTA", "CGAU")
    global w
    w = f'{DNA.translate(trans)}'.strip() #the RNA strand
    return f'The Transcription is: \n 5- {w} -3\n'   
def GC_Content (DNA): 
 a= DNA.count('C') #count C nucleotide in DNA strand
 b= DNA.count('G') #count G nucleotide in DNA strand
 c = ((a+b) / len(DNA)) * 100 #GC contnent in the DNA strand
 return(f"- The GC Content is {c:.1f}%\n")

def protein_translation(RNA):
    "convert RNA sequence to Protein."
    my_codon_dic={ 'UAA':'Phe','UUC':'Phe','UUA':'Leu','UUG':'Leu','CUU':'Leu','CUC':'Leu','CUA':'Leu','CUG':'Leu','AUU':'Ile',
'AUC':'Ile','AUA':'Ile','AUG':'Met',"GUU":'Val','GUC':'Val','GUA':'Val','GUG':'Val','UCU':'Ser',"UCC":'Ser',
'UCA':'Ser','UCG':'Ser','CCU':'Pro','CCC':'Pro',"CCA":'Pro','CCG':'Pro','ACU':'Thr','ACC':'Thr',
'ACA':'Thr','ACG':'Thr',"GCU":'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala','UAU':'Tyr',"UAC":'Tyr','CAU':'His','CAC':'His',
'CAA':'Gln','CAG':'Gln',"AAU":'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys','GAU':'Asp',"GAC":'Asp','GAA':'Glu','GAG':'Glu',
'GGU':'Gly','GGC':'Gly',"GGA":'Gly','GGG':'Gly','AGG':'Arg','AGA':'Arg','AGC':'Ser',"AGU":'Ser','CGG':'Arg','CGA':'Arg',
'CGC':'Arg','CGU':'Arg','UGG':'Trp','UGC':'Cys','UGU':'Cys',"UAA": "STOP","UAG": "STOP","UGA": "STOP",
} #Amino Acids Dictionary 
    protein_strand = [my_codon_dic.get(RNA[i : i + 3],RNA[i : i + 3]) for i in range(0, len(RNA),3)]
    return (f'The Translated Protein A.As is: {protein_strand}\n')
    # use list to put the product of A.As in it
    # we use get method with 2 parameters,the first to give the A.As from the dictionary and the defualt to return the RNA nucleotides if not found in the dictionary.
My_File = open(r'C:\Users\af\Downloads\Documents\python\cmder\Bio Informatics.txt','a')
My_File.write(x)
My_File.write(complementry(DNA))
My_File.write(transcription(strand))
My_File.write(GC_Content(DNA))
My_File.write(f"- The no.of 'A' bases is {DNA.count('A')} and the perecent is {DNA.count('A')/len(DNA)*100:.1f}%\n")
My_File.write(f"- The no.of 'T' bases is {DNA.count('T')} and the perecent is {DNA.count('T')/len(DNA)*100:.1f}%\n")
My_File.write(protein_translation(w))
############################################### Another_Method ################################################################
DNA = input('Please put your sequence:').lower()

def ds(sequence):
    complementary_DNA = sequence.replace("a", "T").replace("c","G").replace("t","A").replace("g","C")
    # complementary_DNA = complementary_DNA.upper()
    print("\nThe double stranded DNA is: \n5' %s 3'" %(sequence.upper()))
    print("3' %s 5'" %(complementary_DNA))

def trnscribed_DNA(sequence):
    RNA = sequence.replace("t", "U")
    print("\nThe RNA sequence is:\n5' %s 3'" %(RNA.upper()))


def translate_dna(sequence):
    """
    translate_dna is a function that can translate any DNA sequence into protein sequence
    and can translate any gap or N nucleotide into "X" residue assuming that these are
    real 1-nucleotide deletions from an in-frame (frame +1) amino acid coding alignment.
    """

    codontable = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
              'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
              'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
              'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
              'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
              'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
              'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
              'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
              'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
              'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
              'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
              'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
              'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
              'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
              'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
              'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
              '---': '-'} # to translate any codon containing gaps only into "-" residue

    protein = ""

    for nt in range(0, len(sequence), 3):
        if sequence[nt:nt + 3] in codontable:
            residue = codontable[sequence[nt:nt + 3]]
        else:
            residue = "-"  #to translate any codon containing gaps or unknown bases like N bases into "X" residue

        protein += residue

    return protein

ds(DNA)

trnscribed_DNA(DNA)

print("\nThe translated protein is : %s" %(translate_dna(DNA.upper())))

print("\nThe base percentage is:\nC: (%i) %.2f%%\nG: (%i) %.2f%%\nA: (%i) %.2f%%\nT: (%i) %.2f%%\n" %(DNA.count('c'), DNA.count('c') * 100 / len(DNA), DNA.count('g'), DNA.count('g') * 100 / len(DNA), DNA.count('a'), DNA.count('a') * 100 / len(DNA), DNA.count('t'), DNA.count('t') * 100 / len(DNA)))DNA = input('Please put your sequence:').lower()

def ds(sequence):
    complementary_DNA = sequence.replace("a", "T").replace("c","G").replace("t","A").replace("g","C")
    # complementary_DNA = complementary_DNA.upper()
    print("\nThe double stranded DNA is: \n5' %s 3'" %(sequence.upper()))
    print("3' %s 5'" %(complementary_DNA))

def trnscribed_DNA(sequence):
    RNA = sequence.replace("t", "U")
    print("\nThe RNA sequence is:\n5' %s 3'" %(RNA.upper()))

def translate_dna(sequence):
    """
    translate_dna is a function that can translate any DNA sequence into protein sequence
    and can translate any gap or N nucleotide into "X" residue assuming that these are
    real 1-nucleotide deletions from an in-frame (frame +1) amino acid coding alignment.
    """

    codontable = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
              'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
              'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
              'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
              'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
              'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
              'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
              'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
              'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
              'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
              'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
              'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
              'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
              'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
              'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
              'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
              '---': '-'} # to translate any codon containing gaps only into "-" residue

    protein = ""

    for nt in range(0, len(sequence), 3):
        if sequence[nt:nt + 3] in codontable:
            residue = codontable[sequence[nt:nt + 3]]
        else:
            residue = "-"  #to translate any codon containing gaps or unknown bases like N bases into "X" residue

        protein += residue

    return Protein

ds(DNA)

trnscribed_DNA(DNA)

print("\nThe translated protein is : %s" %(translate_dna(DNA.upper())))

print("\nThe base percentage is:\nC: (%i) %.2f%%\nG: (%i) %.2f%%\nA: (%i) %.2f%%\nT: (%i) %.2f%%\n" %(DNA.count('c'), DNA.count('c') * 100 / len(DNA), DNA.count('g'), DNA.count('g') * 100 / len(DNA), DNA.count('a'), DNA.count('a') * 100 / len(DNA), DNA.count('t'), DNA.count('t') * 100 / len(DNA)))