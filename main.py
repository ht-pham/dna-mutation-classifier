import pandas as pd
import numpy as np
# import Bio
# from Bio import SeqIO
import random
import re

from classes.dna import DNA
from classes.lambda_dna import LambdaGenome


if __name__ == "__main__":
    # Parent class
    # Example: a randomly generated DNA sequence of random length
    random_dna = DNA()
    random_len = random.randint(0,10000)
    random_seq = random_dna.generate(random_len)
    print(f"Generated random DNA sequence: {random_seq}")
    
    # Protein synthesis: Transcribe DNA to mRNA, then translate to protein
    original_protein = random_dna.synthesizeProtein(random_seq)
    
    # DNA Mutation: computationally mutate DNA from the generated DNA
    # Case 1: Random case
    substitutions,original_seq, mutated_seq = random_dna.alterSNP('random',random_seq)
    print(f"SNP substitution types: {substitutions}")
    aa_mutations = random_dna.findMutations(original_seq , mutated_seq)
    print(aa_mutations)
    print()
    # Protein synthesis on the mutated DNA:
    mutated_protein = random_dna.synthesizeProtein(mutated_seq)

    # Inspect Protein Change
    print(f"Original protein: {original_protein}")
    print(f"Mutated protein: {mutated_protein}")
    protein_changed = 'no' if original_protein == mutated_protein else 'yes'
    print(f"Was protein changed?: {protein_changed} ")
    print('-'*100)

    # Case 2: All case
    all_case_seq = random_dna.generate(random.randint(0,10000))
    print(f"Generated random DNA sequence: {all_case_seq}")
    
    # Protein synthesis: Transcribe DNA to mRNA, then translate to protein
    original_protein = random_dna.synthesizeProtein(all_case_seq)
    
    # DNA Mutation: computationally mutate DNA from the generated DNA
    substitutions,original_seq, mutated_seq = random_dna.alterSNP('all',all_case_seq)
    # Protein synthesis on the mutated DNA:
    for mutated_gene in mutated_seq:
        print("Original ",end="")
        random_dna.display('dna',all_case_seq)
        print("Mutated ",end="")
        mutated_protein = random_dna.synthesizeProtein(mutated_gene)
        
        # Inspect Protein Change
        print(f"Original protein: {original_protein}")
        print(f"Mutated protein: {mutated_protein}")
        protein_changed = 'no' if original_protein == mutated_protein else 'yes'
        print(f"Was protein changed?: {protein_changed} ")
        print('-'*100)

    
    
    # Test set 1: Immunitydeficiency45 / Gene ID: 3455
    header,ifnar2_gene = random_dna.load_genome("data/IFNAR2_datasets/ifnar2_gene.fna")
    
    headers = re.split(r'[>, ]',header)
    #print(headers)
    ifnar2_headers = {"Access ID":headers[1],"Name":headers[2],"Description":" ".join(headers[3:])}
    for k,v in ifnar2_headers.items():
        print(f"{k}: {v}")

    #print(ifnar2_gene[35727:35818]) # second header
    ifnar2_gene = "".join(ifnar2_gene[:35727]+ifnar2_gene[35818:])
    print(f"First 102-nucleotide base of {ifnar2_headers['Name']} genome")
    random_dna.display('dna',sequence=ifnar2_gene[:102])

    chromosome_21 = random_dna.countAll('dna',ifnar2_gene)
    first = int(input("Enter the first codon position of the subsequence: "))
    last = int(input("Enter the last codon position of the subsequence: "))
    subsequence = random_dna.getSubsequence(first,last)
    random_dna.display('dna',sequence=subsequence)
    print("_"*100)
    print()
    # Explore the first coding region and its rna
    print("DNA, RNA, and Protein sequence of the first coding region of IFNAR2 gene (Chromosome 21):")
    # random_dna.display('dna',ifnar2_gene[585:624])
    # print("Its corresponding RNA:")
    # rna_ifnar2_first = random_dna.transcribe(ifnar2_gene[585:624])
    # random_dna.display('rna',rna_ifnar2_first)
    # first_protein = random_dna.translate(rna_ifnar2_first)
    # print(f"Protein: {first_protein}")
    ifnar2_first_protein = random_dna.synthesizeProtein(ifnar2_gene[585:624])
    
    rna_ifnar2 = random_dna.transcribe(ifnar2_gene)
    #random_dna.display('RNA',rna_ifnar2)

    
    print("Look up codons on RNA sequence of IFNAR2 gene (Chromosome 21):")
    chromosome_21 = random_dna.countAll('rna',rna_ifnar2)
    prompt = input("Enter codon that you want to look up: (e.g. UAA, AUG, etc.): ")
    print(random_dna.getPositionList(prompt))
    
    first = int(input("Enter the first codon position of the subsequence: "))
    last = int(input("Enter the last codon position of the subsequence: "))
    subsequence = random_dna.getSubsequence(first,last)
    random_dna.display('rna',sequence=subsequence)
    
    print("_"*100)
    print()
    print("Full RNA of IFNAR2 genome")
    random_dna.display('RNA',rna_ifnar2)
    print("_"*100)
    print()
    
    # Test set 2: Lambda phage
    # Child class
    lambda_dna = LambdaGenome("Enterobacteria phage lambda",None) # create object for Lambda Genome class
    
    # load .fasta file into a tuple (metadata, dna_sequence)
    header,lambda_genome = lambda_dna.load_genome()
    
    # Get general info of the data from the .fasta file
    headers = re.split(r'[>, ]',header)
    lambda_headers = {"Access ID":headers[1],"Name":" ".join(headers[2:5]),"Description":" ".join(headers[5:])}
    for k,v in lambda_headers.items():
        print(f"{k}: {v}")
    #lambda_dna.printStatistics('dna',lambda_genome)
    
    # Display 34 first codons
    print(f"First 34 codons",end=" ")
    lambda_dna.display('dna',lambda_genome[:102])
    
    # Display first coding region of mRNA
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
    subsequence=lambda_dna.getSubsequence(first,last)
    lambda_dna.display('dna',sequence=subsequence)

    # Look up 3: Find all positions of a specific codon on the DNA
    mrna = lambda_dna.transcribe(lambda_genome)
    codons = lambda_dna.countAll('rna',mrna)
    prompt = input("Enter codon that you want to look up: (e.g. UAA, AUG, etc.): ")
    print(lambda_dna.getPositionList(prompt))
    
    
    
