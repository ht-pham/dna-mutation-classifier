import pandas as pd
import numpy as np
# import Bio
# from Bio import SeqIO
import random

from classes.dna import DNA
from classes.lambda_dna import LambdaGenome


if __name__ == "__main__":
    # Parent class
    # Example: a randomly generated DNA sequence of random length
    random_dna = DNA()
    random_len = random.randint(0,10000)
    random_seq = random_dna.generate(random_len)
    print(f"Generated random DNA sequence: {random_seq}")
    
    # Break the randomly generated DNA sequence down to codons and their positions in the sequence
    #codons = random_dna.decipher(random_seq)
    #random_dna.printStatistics('dna',random_seq)
   
    # Protein synthesize Transcribe to mRNA, then translate to protein
    original_protein = random_dna.synthesizeProtein(random_seq)
    
    # Syntheticize the generated DNA (i.e. computationally generate mutated DNA)
    substitutions,_ , mutated_seq = random_dna.alterSNP(random_seq)
    #print(f"Mutated DNA sequence: {mutated_seq}")
    print(f"SNP substitution types: {substitutions}")
    #random_dna.printStatistics('dna',mutated_seq)
    mutated_protein = random_dna.synthesizeProtein(mutated_seq)
    
    print(f"Original protein: {original_protein}")
    print(f"Mutated protein: {mutated_protein}")
    protein_changed = 'no' if original_protein == mutated_protein else 'yes'
    print(f"Was protein changed?: {protein_changed} ")
    print('-'*100)
    
    # Child class
    lambda_dna = LambdaGenome("Enterobacteria phage lambda",None) # create object for Lambda Genome class
    # load .fasta file into a tuple (metadata, dna_sequence)
    header,lambda_genome = lambda_dna.load_genome()
    # Get general info of the data
    import re
    headers = re.split(r'[>, ]',header)
    lambda_headers = {"Access ID":headers[1],"Name":" ".join(headers[2:5]),"Description":" ".join(headers[5:])}
    for k,v in lambda_headers.items():
        print(f"{k}: {v}")
    # left_cos_site = lambda_genome[0:11]
    # print(f"Left cos site",end=" ")
    # lambda_dna.display('dna',left_cos_site)
    print(f"First 34 codons",end=" ")
    lambda_dna.display('dna',lambda_genome[:102])
    #lambda_dna.printStatistics('dna',lambda_genome)
    print(f"mRNA of Lambda Phage (only first coding region)",end=" ")
    lambda_dna.display('rna',lambda_dna.transcribe(lambda_genome)[:51])
    #lambda_dna.printStatistics('rna',lambda_dna.transcribe(lambda_genome))
    
    codons = lambda_dna.countAll('dna',lambda_genome)
    # Look up 1: Find codon at specific order
    prompt = int(input("Enter the order (1-16168) of the codon: "))
    print(lambda_dna.getAminoAcidAt(prompt))

    # Look up 2: Look up subsequence based on the order of codon
    first = int(input("Enter the order (1-16168) of the first codon: "))
    last = int(input("Enter the order (1-16168) of the last codon: "))
    print(lambda_dna.getSubsequence(first,last))

    # Look up 3: Find all positions of a specific codon on the DNA
    mrna = lambda_dna.transcribe(lambda_genome)
    codons = lambda_dna.countAll('rna',mrna)
    prompt = input("Enter codon that you want to look up: (e.g. UAA, AUG, etc.): ")
    print(lambda_dna.getPositionList(prompt))
    
    
    