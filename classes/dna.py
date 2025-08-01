class DNA:
    
    def __init__(self):
        # 4 bases/nucleotides
        self.nucleotides = ['A','T','C','G']
        # Mapping rules for each nucleotides
        self.genetic_codes = dict(zip(self.nucleotides,[['A','T'],['T','A'],['C','G'],['G','C']]))
        # Amino Acid (Type/Name) Chart
        self.rna_amino_acids_chart = {'phe':['UUU','UUC'],'leu':['UUA','UUG','CUU','CUC','CUA','CUG'],
                       'ser':['UCU','UCC','UCA','UCG','AGU','AGC'],'tyr':['UAU','UAC'],
                       'cys':['UGU','UGC'],'try':'UGG','pro':['CCU','CCC','CCA','CCG'],
                       'stop':['UAA','UAG','UGA'],'his':['CAU','CAC'],'gln':['CAA','CAG'],
                       'arg':['CGU','CGC','CGA','CGG','AGA','AGG'],'lys':['AAA','AAG'],
                       'asp':['AAU','AAC','GAU','GAC'],'thr':['ACU','ACC','ACA','ACG'],
                       'ile':['AUU','AUC','AUA'],'met':'AUG','glu':['GAA','GAG'],
                       'val':['GUU','GUC','GUA','GUG'],'ala':['GCU','GCC','GCA','GCG'],
                       'gly':['GGU','GGC','GGA','GGG']
                       } 
        self.dna_amino_acids_chart = {'phe':['TTT','TTC'],'leu':['TTA','TTG','CTT','CTC','CTA','CTG'],
                       'ser':['TCT','TCC','TCA','TCG','AGT','AGC'],'tyr':['TAT','TAC'],
                       'cys':['TGT','TGC'],'try':'TGG','pro':['CCT','CCC','CCA','CCG'],
                       'stop':['TAA','TAG','TGA'],'his':['CAT','CAC'],'gln':['CAA','CAG'],
                       'arg':['CGT','CGC','CGA','CGG','AGA','AGG'],'lys':['AAA','AAG'],
                       'asp':['AAT','AAC','GAT','GAC'],'thr':['ACT','ACC','ACA','ACG'],
                       'ile':['ATT','ATC','ATA'],'met':'ATG','glu':['GAA','GAG'],
                       'val':['GTT','GTC','GTA','GTG'],'ala':['GCT','GCC','GCA','GCG'],
                       'gly':['GGT','GGC','GGA','GGG']
                       } 
        # Basic codons
        self.start = 'ATG'
        self.stop = ['TAG','TGA','TAA']

        self.sequence = self.start

        self.codons = []
        self.amino_acids = set()
        
    def __str__(self):
        return f"{self.sequence} has length of {len(self.sequence)}"

    def generate(self,num_nucleotides):
        '''
            Input: Length of the random DNA sequence
            Output: Randomly generated DNA sequence with start codon as 'ATG'
        '''
        import random
        self.sequence = self.start
        choice = ""
        for i in range(num_nucleotides):
            last_nucleotide = self.sequence[-1:]
            if i%3 == 0:
                choice = random.choice(self.nucleotides)
            else:
                choice = random.choice(self.genetic_codes[last_nucleotide])
            self.sequence = self.sequence + choice
                 
            # if STOP codon has been found
            last_codon = self.sequence[-3:] == 'TAG' or self.sequence[-3:] == 'TGA' or self.sequence[-3:] == 'TAA'
            #print(last_codon)
            if last_codon:
                break
        
        return self.sequence
    
    def load_genome(self,file_path):
        '''
            Input: Fasta file's directory
            Output: A tuple of .fasta header and the complete DNA sequence
        '''
        genome = []
        with open(file_path,'r') as f:
            for line in f:
                current_line = line.strip()
                genome.append(current_line)

        # Separate file header from the loaded strings
        header = genome[0] 
        # Merge the rest of the list into 1 string to make a complete DNA sequence
        self.sequence = "".join(genome[1:])

        return header,self.sequence

    def getSequence(self):
        return self.sequence
    
    def decipher(self,sequence):
        '''
            Input: a DNA/RNA sequence as string type
            Output: a dictionary of 3-base codons in the DNA sequence
        '''
        #self.sequence = dna_sequence
        codons = []
        for i in range(0,len(sequence),3):
            codons.append(sequence[i:i+3])
            #print(f"{dna_sequence[i:i+3]}",end=" ")
        self.codons = dict(zip(range(1,len(codons)+1),codons))
        return self.codons
    
    def setPositionList(self,codons_dict):
        '''
            Input: a dictionary of positions and their corresponding codon
                    e.g. {1:'ATG',2:'GGT',3:'ATG',4:'TAG'}
            Output: a reversed dictionary of each codon and their position
                    e.g. {'ATG':[1,3],'GGT':[2],'TAG':[4]}
        '''
        self.codons = codons_dict
        # Convert the list of 3-base codons to a set of unique values
        self.amino_acids = set(self.codons.values())

        position_list = {codon: [] for codon in self.amino_acids}
        for k,v in self.codons.items():
            position_list[v].append(k) # collect positions of each codon

        return position_list
    
    def getPositionList(self,codon):
        '''
            Input: a specific codon e.g. 'GCG'
            Output: the list of all positions/orders of that codon appearing in the DNA
        '''
        try:
            #codons = self.getSequence()
            positions = self.setPositionList(self.codons)
            return positions[codon.upper()]
        
        except KeyError:
            print("Error occured: Wrong typo or Codon does not exist in the genome.")
            return []

    def getSubsequence(self,first,last):
        '''
            Input: Order of the first codon and the last codon
            Output: the sub-sequence of DNA from the first to the last
        '''
        subseq = []
        for i in range(first,last+1):
            subseq.append(self.codons[i])
    
        return "".join(subseq)

    def getAminoAcidAt(self,codon_index):
        '''
            Input: a specific order (of a codon) in DNA sequence
            Output: the codon on that order
        '''
        return self.codons.get(codon_index)
    
    def getAminoAcidName(self,nucleic_acid,codons):
        '''
            Input: a set of codons in the DNA sequence
            Output: A list of Amino Acids' Names based on the Amino Acid's Name Chart
        '''
        table = self.dna_amino_acids_chart if nucleic_acid.upper() == 'dna'.upper() else self.rna_amino_acids_chart
        amino_names = dict(zip(codons,['']*len(codons)))
        
        # Look up amino acids' name in the chart
        for k,v in table.items():
            # for each codon in the DNA sequence
            for codon in codons:
                if str(codon) in v: 
                    amino_names[codon] = k
                elif len(str(codon))<3:
                    amino_names[codon] = 'Incomplete'
        
        return amino_names

    def countAll(self,nucleic_acid,sequence):
        '''
            Input: a DNA/mRNA sequence as string type
            Output: The corresponding amino acids and their occurrence in the DNA sequence
        '''
        self.codons = self.decipher(sequence)
        # Get positions of each codon
        position_list = self.setPositionList(self.codons)
        # Convert the list of 3-base codons to a set of unique values
        self.amino_acids = set(self.codons.values()) 
        amino_names = self.getAminoAcidName(nucleic_acid,self.amino_acids)
        
        stats = {codon: 0 for codon in self.amino_acids}
        for codon in self.codons.values():
            if codon in stats.keys():
                stats[codon]+=1

        return amino_names, stats, position_list
    
    def display(self,nucleic_acid,sequence):
        self.codons = self.decipher(sequence)
        aa = self.getAminoAcidName(nucleic_acid,set(self.codons.values()))
        
        # First start codon flag
        first_start_codon = False
        first_stop_codon = False

        BACKGROUND_BRIGHT_RED = '\033[101m'
        BACKGROUND_BRIGHT_GREEN = '\033[102m'
        BACKGROUND_BRIGHT_WHITE = '\033[107m'
        RED ='\033[31m'
        BOLD_RED = '\u001b[91;1m'
        GREEN = '\033[32m'
        BOLD_GREEN = '\u001b[92;1m'
        RESET = '\033[0m'
        print("Codon Sequence:")
        for i in range(1,len(self.codons.values())+1):
            if i == len(self.codons.values()):
                print(f"[{i}] {BOLD_RED}{self.codons[i]}{RESET}")
            else:
                if (self.codons[i] in self.stop) or (self.codons[i] in ['UAG','UGA','UAA']):
                    if first_start_codon == True and first_stop_codon == False: # First stop codon
                        first_stop_codon = True   
                        print(f"[{i}] {BOLD_RED}{BACKGROUND_BRIGHT_RED}{self.codons[i]}{RESET}",end=" -> ")
                    else:
                        print(f"[{i}] {RED}{self.codons[i]}{RESET}",end=" -> ")
                elif (self.codons[i] == 'ATG') or (self.codons[i] == 'AUG'):
                    if first_start_codon == False:
                        first_start_codon = True
                        print(f"[{i}] {BOLD_GREEN}{BACKGROUND_BRIGHT_GREEN}{self.codons[i]}{RESET}", end=" -> ")
                    else:
                        print(f"[{i}] {GREEN}{self.codons[i]}{RESET}",end=" -> ")
                else:
                    print(f"[{i}] {self.codons[i]}",end=" -> ")
        
        # First start codon flag
        first_start_codon = False
        first_stop_codon = False
        print("Amino Acids:")
        for pos,codon in self.codons.items():
            if pos == len(self.codons.values()):
                print(f"[{pos}] {RED}{aa[codon]}{RESET}")
            else:
                if aa[codon] == 'stop':
                    if first_stop_codon == False:
                        first_stop_codon = True   
                        print(f"[{i}] {BOLD_RED}{BACKGROUND_BRIGHT_RED}{aa[codon]}{RESET}",end=" -> ")
                    else:
                        print(f"[{pos}] {RED}{aa[codon]}{RESET}",end=" -> ")
                elif aa[codon] == 'met':
                    if first_start_codon == False:
                        first_start_codon = True
                        print(f"[{i}] {BOLD_GREEN}{BACKGROUND_BRIGHT_GREEN}{aa[codon]}{RESET}", end=" -> ")
                    else:
                        print(f"[{i}] {GREEN}{aa[codon]}{RESET}",end=" -> ")
                else:
                    print(f"[{pos}] {aa[codon]}",end=" -> ")

    def printStatistics(self,nucleic_acid,sequence):
        print("Statistics:")
        names, stats, position_list = self.countAll(nucleic_acid,sequence) 
        print("-"*62)
        print("| Codon | Amino Acid | Occurrence |          Order(s)         |")
        print("-"*62)
        for k,v in stats.items():
            if len(position_list.get(k))<10:
                pos_str = ", ".join(str(p) for p in position_list.get(k))
                print(f"|  {str(k):^3}  | {names.get(k):^10} |   {v:^7}  |     {pos_str:^12}          |")
            else:
                first = position_list.get(k)[0]
                last = position_list.get(k)[-1]
                print(f"|  {str(k):^3}  | {names.get(k):^10} |   {v:^7}  | First: {first:<3} & Last: {last:>3}  |")
        print("-"*62)
    
    def transcribe(self,sequence):
        '''
            Input: a DNA sequence
            Output: a list of mRNAs
            Approach: Using Open Reading Frame (ORF) analysis
        '''
        # First, break down the sequence into codons

        codons = self.decipher(sequence.upper()) 
        
        # Second, define coding segments in the DNA sequence
        start_pos = self.getPositionList('ATG')
        #print(start_pos) # check if any
        if start_pos == []:
            print("Only non-coding region given: no start codon found")
            print("Abort DNA transcription.")
            return sequence
        
        try:
            # Original Approach
            #---------------------------------
            # stop_pos = []
            # for pos,codon in codons.items():
            #     if codon in self.stop:
            #         stop_pos.append(pos)
            #---------------------------------
            stop_pos = []
            
            # 2nd approach
            # Choose the codon after the first start codon
            pointer = start_pos[0]+1
            index = 0
            # Run from that codon to the last codon of the sequence
            while pointer < len(codons.values()):
                codon = codons[pointer]
                if pointer in start_pos:
                    pointer = pointer+1
                    continue
                if codon in self.stop and pointer not in stop_pos:
                    stop_pos.append(pointer)
                    # Move to the next start codon
                    index = index + 1
                    if index >= len(start_pos):
                        break
                    else: 
                        if start_pos[index] < pointer:
                            # Move to the next start codon
                            index = index + 1
                            pointer = start_pos[index]    
                pointer = pointer+1
            
            # Test duplicates
            # print(len(stop_pos)) 
            # print(len(set(stop_pos)))

            # if no stop codon in the DNA sequence
            if stop_pos == []:
                stop_pos = [len(codons.values())]
        
            #print(stop_pos)

            # define coding regions
            # Original Approach: Impractical
            # coding_regions = dict(zip(start_pos,stop_pos))
            # 2nd approach: 
            # (1) Combine both lists into one list, 
            # (2.1) remove the 'ATG' positions in between 'ATG' and stop codons
            # (2.2) remove the stop codons' position following after stop codons (i.e. no start codon in between them)
            # (3) Pair two adjacent elements to one dictionary
            start_stop = sorted(start_pos+stop_pos)
            #print(start_stop)
            
            is_coding = False # Flag to remove repeated start and stop codons
            
            for index,pos in enumerate(start_stop):
                # if the previous codon is a stop codon and the current codon is a start
                # (i.e. in non-coding region => in coding region)
                if pos in start_pos and is_coding==False:
                    # then change the flag
                    is_coding = True
                # if the previous codon is a start codon and the current codon is a start 
                # (i.e. already in coding region => mark as 0)   
                elif pos in start_pos and is_coding==True:
                    # then pop the current codon out of the list
                    start_stop[index] = 0
                # if the current codon is a stop codon in coding region, then change the flag to False
                # (i.e. in coding region => in non-coding region)
                elif pos in stop_pos and is_coding==True:
                    is_coding = False
                # if the previous codon is a stop codon and the current codon is a stop codon
                # (i.e. already in non-coding region => mark as 0)
                elif pos in stop_pos and is_coding==False:
                    start_stop[index] = 0
            
            start_stop = [x for x in start_stop if x != 0] # remove the repeated codons
            #print(len(start_stop))
            #print(start_stop)
            coding_regions = dict(zip(start_stop[::2],start_stop[1::2]))

            #print(coding_regions)

            # Third, retrieve coding DNA sequence from start and stop positions
            pre_mRNAs = []
            for start,stop in coding_regions.items():
                pre_mRNA = self.getSubsequence(start,stop)
                pre_mRNA = pre_mRNA.replace("T","U")
                pre_mRNAs.append(pre_mRNA)
            
            mRNAs = "".join(pre_mRNAs) # this is biologically incorrect because we need to remove non-coding regions within each pre_mRNA before generating mRNA
        
            return mRNAs
        
        except KeyError:
            print("Start Loss Mutation Occured")
            return sequence


    def translate(self,mRNA):
        # When start loss mutation occured, it likely already failed to generate mRNA
        # return the DNA sequence without translation
        if 'T' in mRNA:
            print("Only non-coding region: Failed to create protein")
            return mRNA
        
        codons = self.decipher(mRNA)
        amino_acids = self.getAminoAcidName('rna',codons.values())
        for codon,aa in amino_acids.items():
            if aa == 'stop':
                amino_acids[codon] = "*"
            elif aa == 'Incomplete':
                amino_acids[codon] = ''
            else:
                #amino_acids[codon] = amino_acids[codon]+"-"
                amino_acids[codon]= aa[0].upper()+aa[1:]
        
        aa = []
        for codon in codons.values():
            aa.append(amino_acids[codon])

        proteins ="".join(aa)
        return proteins
    
    def synthesizeProtein(self,dna_sequence):
        #print(f"DNA: {dna_sequence}")
        print("DNA",end=" ")
        self.display('dna',dna_sequence)

        #transcribe DNA to mRNA
        mRNA = self.transcribe(dna_sequence)
        if 'T' in mRNA:
            #this means, start loss mutation occured. Failed to create protein
            return "Incomplete"
        
        print("mRNA",end=" ")
        self.display('rna',mRNA)
        #print(f"Coresponding mRNA: {mRNA}")
        #self.printStatistics('rna',mRNA)
        
        #translate mRNA to Protein
        protein = self.translate(mRNA)
        print(f"Protein: {protein}")
        print("_"*100)
        print()

        return protein
    
    
    def alterSNP(self,type_of_SNP='random',original_sequence='ATGTAG'):
        '''
            Description: A synthetic Single Nucleotide Polymorphism function 
            Input: a DNA sequence
            Output: the original sequence and its mutated sequence
        '''
        def substitute(genome,snp_location,new_nucleotide):
            # Classify bases based on its chemical structures
            heterocyclic_bases = { 'purine':{'A','G'},'pyrimidine':{'T','C'}}
            substitution_types = []

            if new_nucleotide != genome[snp_location]:
                bases = set([genome[snp_location],new_nucleotide])
                
                # highlight the base that gets changed
                RED ='\033[91m'
                RESET = '\033[0m'
                genome_string_before = "".join(genome[:snp_location])
                genome_string_after = "".join(genome[(snp_location+1):])
                print(f"Original DNA: {genome_string_before}[{RED}{genome[snp_location]}{RESET}]{genome_string_after}")
                #print(f'Mutation at {snp_location+1}: {genome[snp_location]} -> {new_nucleotide}')
                
                # Replace the selected single base with the new base
                genome[snp_location] = new_nucleotide

                # highlight the base that got changed
                genome_string_before = "".join(genome[:snp_location])
                genome_string_after = "".join(genome[(snp_location+1):])
                print(f"Mutated DNA: {genome_string_before}[{RED}{genome[snp_location]}{RESET}]{genome_string_after}")
                print()
                ## then, label SNP substitution
                #### Transition: purine <-> purine, pyrimidine <-> pyrimidine
                #### Transversion: purine <-> pyrimidine 
                same_heterocyclic_base = bases in heterocyclic_bases.values()
                if same_heterocyclic_base:
                    substitution_types.append('transition')
                else:
                    substitution_types.append('transversion')
            
            mutated_genome = "".join(genome)
            
            return substitution_types, mutated_genome
        
        import math

        # Convert the sequence for the alternation step
        genome = list(original_sequence.upper()) # e.g. 'ATGTAG' => ['A', 'T', 'G', 'T', 'A', 'G']

        # Decide number of mutation randomly (it cannot exceed length of DNA sequence)
        import random
        times = random.randint(1,len(genome))
        #print(f"Potential number of mutations: {times}")

        snp_location = 0 # start from the first nucleotide
        substitution_types = []
        
        while times>0:
            # First, select a random single base from the sequence (Starting point at index 0)
            snp_location = snp_location + random.randint(1,len(original_sequence)-1) # first loop: 0 + random int
            # if the random location surpasses the length of the DNA sequence, then break the loop
            if snp_location >= len(original_sequence):
                break
            codon_pos = math.floor(snp_location/3)+1
            # else move forward with single nucleotide mutation
            # Second, pick a random single nucleotide if type_of_SNP is 'random'
            if type_of_SNP=='random':
                random_nucleotide = random.choice(['A','T','C','G'])
                substitution_types,mutated_genome = substitute(genome,snp_location,random_nucleotide)
            else:
                all_SNPs = {}
                all_SNPs[snp_location]=[]
                mutated_genomes = []
                substitution_types = []
                for i,nucleotide in enumerate(self.nucleotides):
                    substitution_type,mutated_genome = substitute(list(original_sequence.upper()),snp_location,nucleotide)
                    aa_mutations = self.findMutations(original_sequence=original_sequence,mutated_genome=mutated_genome)
                    if substitution_type != []:
                        all_SNPs[snp_location].append([substitution_type,mutated_genome,aa_mutations[codon_pos]])
                        mutated_genomes.append([snp_location,i+1,mutated_genome])
                        substitution_types.append([snp_location,i+1,substitution_type])
                    
                mutated_genome = [x[2] for x in mutated_genomes]
                #print(mutated_genome)
                #print(substitution_types)
                
                # Print out purpose only
                for num,case in enumerate(all_SNPs[snp_location]):
                    print(f"Case {num+1}:")
                    print(f"- Mutated gene: {case[1]}")
                    print(f"- SNP Type: {case[0][0]}")
                    print(f"- [{codon_pos}] {case[2]['codon'][0]} -> {case[2]['codon'][1]}: {case[2]['amino'][0]} -> {case[2]['amino'][1]}")
                    
            # Reduce count
            times = times-1

        
        return substitution_types, original_sequence, mutated_genome

    def findMutations(self,original_sequence, mutated_genome):
        '''
            Description: find the codon position(s) and their original and their new codons and amino acids
            Input: old DNA sequence, new DNA sequence
            Output: a dictionary of codon positions and their changes.
                {pos_mutation:{codon:[original codon,new codon],aa:{original aa, new aa}}}
        '''
        original_codons = self.decipher(original_sequence)
        self.codons = {} # reset self.codons values
        mutated_codons = self.decipher(mutated_genome)
        
        pos_mutations = []
        for codon_pos in original_codons.keys():
            if original_codons[codon_pos] != mutated_codons[codon_pos]:
                pos_mutations.append(codon_pos)
        # print(pos_mutations)
        
        # dict = {pos_mutation:{codon:[original codon,new codon],aa:{original aa, new aa}}}
        mutations = {}
        for codon_pos in pos_mutations:
            old = original_codons[codon_pos]
            new = mutated_codons[codon_pos]
            aa = self.getAminoAcidName('dna',[old,new])
            mutations[codon_pos]={'codon':[old,new],'amino':list(aa.values())}
    
        #print(mutations)
        return mutations
        
        
    # def labelMutation(self,mutation):
    #     mutation_type = ["synonymous","missense","nonsense"]