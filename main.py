import pandas as pd
import numpy as np
# import Bio
# from Bio import SeqIO
import random

from classes.dna import DNA
from classes.lambda_dna import LambdaGenome


if __name__ == "__main__":
    # Parent class
    # Example: a randomly generated DNA sequence
    random_dna = DNA()
    random_len = random.randint(0,10000)
    random_seq = random_dna.generate(random_len)
    codons = random_dna.decipher(random_seq)
    print(f"Generated random DNA sequence: {random_seq}")
    random_dna.printStatistics('dna',random_seq)
    # Transcript the random DNA
    print("Transcripting DNA to mRNA...")
    mRNA = random_dna.transcript(random_seq)
    print(f"mRNA: {mRNA}")
    random_dna.printStatistics('rna',random_dna.transcript(random_seq))
    
    # Translate the mRNA to protein
    protein = random_dna.translate(mRNA)
    print(f"Protein: {protein}")
    
    print("_"*100)
    print()
    
    
    # Child class
    lambda_dna = LambdaGenome("Enterobacteria phage lambda",None) # create object for Lambda Genome class
    # gene = lambda_dna.load_from_file()  # load the .fasta file to assign the complete gene to the 'gene' object's sequence  
    
    # load .fasta file into a tuple (metadata, dna_sequence)
    header,lambda_genome = lambda_dna.load_genome()
    # Get general info of the data
    import re
    headers = re.split(r'[>, ]',header)
    lambda_headers = {"Access ID":headers[1],"Name":" ".join(headers[2:5]),"Description":" ".join(headers[5:])}
    for k,v in lambda_headers.items():
        print(f"{k}: {v}")
    left_cos_site = lambda_genome[0:11]
    print(f"Left cos site: {left_cos_site}")
    print(f"Sample: {lambda_genome[:100]}")
    lambda_dna.printStatistics('dna',lambda_genome)
    print(f"mRNA of Lambda Phage: \n{lambda_dna.transcript(lambda_genome)}")
    lambda_dna.printStatistics('rna',lambda_dna.transcript(lambda_genome))

    # Look up 1: Find codon at specific order
    prompt = int(input("Enter the order (1-16168) of the codon: "))
    print(lambda_dna.getAminoAcidAt(prompt))

    # Look up 2: Look up subsequence based on the order of codon
    first = int(input("Enter the order (1-16168) of the first codon: "))
    last = int(input("Enter the order (1-16168) of the last codon: "))
    print(lambda_dna.getSubsequence(first,last))

    first = int(input("Enter the order (1-16168) of the first codon: "))
    last = int(input("Enter the order (1-16168) of the last codon: "))
    print(lambda_dna.getSubsequence(first,last))

    # Look up 3: Find all positions of a specific codon on the DNA
    prompt = input("Enter codon that you want to look up: (e.g. UAA, AUG, etc.): ")
    print(lambda_dna.getPositionList(prompt))
    
    
    