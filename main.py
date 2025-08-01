import pandas as pd
import numpy as np
# import Bio
# from Bio import SeqIO
import random
import re

from classes.dna import DNA
from classes.lambda_dna import LambdaGenome

def run_random_case(dna=DNA(),sequence_len=random.randint(0,10000)):
    # Example: a randomly generated DNA sequence of random length
    random_seq = dna.generate(sequence_len)
    print(f"Generated random DNA sequence: {random_seq}")
    dna.display('dna',random_seq)
    
    # Protein synthesis: Transcribe DNA to mRNA, then translate to protein
    original_protein = dna.synthesizeProtein(random_seq)
    # Case 1: Random case
    substitutions,original_seq, mutated_seq = dna.alterSNP('random',random_seq)
    print(f"SNP substitution types: {substitutions}")
    aa_mutations = dna.findMutations(original_seq , mutated_seq)
    print(aa_mutations)
    print()
    # Protein synthesis on the mutated DNA:
    mutated_protein = dna.synthesizeProtein(mutated_seq)

    # Inspect Protein Change
    print(f"Original protein: {original_protein}")
    print(f"Mutated protein: {mutated_protein}")
    protein_changed = 'no' if original_protein == mutated_protein else 'yes'
    print(f"Was protein changed?: {protein_changed} ")
    print('-'*100)

def run_all_cases(dna=DNA(),sequence_len=random.randint(0,10000)):
    all_case_seq = dna.generate(sequence_len)
    print(f"Generated random DNA sequence: {all_case_seq}")
    
    # Protein synthesis: Transcribe DNA to mRNA, then translate to protein
    original_protein = dna.synthesizeProtein(all_case_seq)
    
    # DNA Mutation: computationally mutate DNA from the generated DNA
    substitutions,original_seq, mutated_seq = dna.alterSNP('all',all_case_seq)
    # Protein synthesis on the mutated DNA:
    for i,mutated_gene in enumerate(mutated_seq):
        print("Original ",end="")
        dna.display('dna',original_seq)
        print("Mutated ",end="")
        mutated_protein = dna.synthesizeProtein(mutated_gene)

        print(f"SNP Type: {substitutions[i][2][0]} at nucleotide #{substitutions[i][0]+1}")
        # Inspect Protein Change
        print(f"Original protein: {original_protein}")
        print(f"Mutated protein: {mutated_protein}")
        protein_changed = 'no' if original_protein == mutated_protein else 'yes'
        print(f"Was protein changed?: {protein_changed} ")
        print('-'*100) 
    
def load_and_run(dna=DNA(),file_path="data/covid.txt"):
    header,sequence = dna.load_genome(file_path)
    headers = re.split(r'[>, ]',header)
    header_df = {"Access ID":headers[1],"Name":headers[2],"Description":" ".join(headers[3:])}
    #for k,v in header_df.items():
    #    print(f"{k}: {v}")

    return header_df,sequence

def covid_mini_program():
    header,covid=load_and_run(DNA(),"data/covid.txt")
    print(f"Access ID: {header['Access ID']}, Name & Description: {header['Name']} {header['Description']}")
    #DNA().display('dna',covid[900:2000])
    rna = DNA().synthesizeProtein(covid[900:2000])

def ifnar2_mini_program():
    # Test set 1: Immunitydeficiency45 / Gene ID: 3455
    dna=DNA()
    header,ifnar2_gene = load_and_run(dna,"data/IFNAR2_datasets/ifnar2_gene.fna")
    print(f"Access ID: {header['Access ID']}, Name & Description:  {header['Name']} {header['Description']}")

    ifnar2_gene = "".join(ifnar2_gene[:35727]+ifnar2_gene[35818:])
    print(f"First 102-nucleotide base of {header['Name']} genome")
    dna.display('dna',sequence=ifnar2_gene[:102])

    chromosome_21 = dna.countAll('dna',ifnar2_gene)
    first = int(input("Enter the first codon position of the subsequence: "))
    last = int(input("Enter the last codon position of the subsequence: "))
    subsequence = dna.getSubsequence(first,last)
    dna.display('dna',sequence=subsequence)
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
    ifnar2_first_protein = dna.synthesizeProtein(ifnar2_gene[585:624])
    
    rna_ifnar2 = dna.transcribe(ifnar2_gene)
    #random_dna.display('RNA',rna_ifnar2)

    
    print("Look up codons on RNA sequence of IFNAR2 gene (Chromosome 21):")
    chromosome_21 = dna.countAll('rna',rna_ifnar2)
    prompt = input("Enter codon that you want to look up: (e.g. UAA, AUG, etc.): ")
    print(dna.getPositionList(prompt))
    
    first = int(input("Enter the first codon position of the subsequence: "))
    last = int(input("Enter the last codon position of the subsequence: "))
    subsequence = dna.getSubsequence(first,last)
    dna.display('rna',sequence=subsequence)
    
    print("_"*100)
    print()
    #print("Full RNA of IFNAR2 genome")
    #dna.display('RNA',rna_ifnar2)
    

def lambda_mini_program():
    # Child class
    lambda_dna = LambdaGenome("Enterobacteria phage lambda",None) # create object for Lambda Genome class
    
    header,lambda_phage = load_and_run(lambda_dna,"data/lambda_genome.fasta")
    print(f"Access ID: {header['Access ID']}, Name & Description:  {header['Name']} {header['Description']}")

    # Display 34 first codons
    print(f"First 34 codons",end=" ")
    lambda_dna.display('dna',lambda_phage[:102])
    
    # Display first coding region of mRNA
    print(f"mRNA of Lambda Phage (only first coding region)",end=" ")
    lambda_dna.display('rna',lambda_dna.transcribe(lambda_phage)[:51])
    #lambda_dna.printStatistics('rna',lambda_dna.transcribe(lambda_genome))
    
    codons = lambda_dna.countAll('dna',lambda_phage)
    # Look up 1: Find codon at specific order
    prompt = int(input("Enter the order (1-16168) of the codon: "))
    print(lambda_dna.getAminoAcidAt(prompt))

    # Look up 2: Look up subsequence based on the order of codon
    first = int(input("Enter the order (1-16168) of the first codon: "))
    last = int(input("Enter the order (1-16168) of the last codon: "))
    subsequence=lambda_dna.getSubsequence(first,last)
    lambda_dna.display('dna',sequence=subsequence)

    # Look up 3: Find all positions of a specific codon on the DNA
    mrna = lambda_dna.transcribe(lambda_phage)
    codons = lambda_dna.countAll('rna',mrna)
    prompt = input("Enter codon that you want to look up: (e.g. UAA, AUG, etc.): ")
    print(lambda_dna.getPositionList(prompt))

if __name__ == "__main__":
    # Parent class
    # Example: a randomly generated DNA sequence of random length
    run_random_case()
    run_all_cases()

    #covid_mini_program()
    ifnar2_mini_program()
    #lambda_mini_program()
    
    
    
