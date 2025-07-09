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
            Output: a list of 3-base codons in the DNA sequence
        '''
        #self.sequence = dna_sequence
        self.codons = []
        for i in range(0,len(sequence),3):
            self.codons.append(sequence[i:i+3])
            #print(f"{dna_sequence[i:i+3]}",end=" ")
        self.codons = dict(zip(range(1,len(self.codons)+1),self.codons))
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
        
        positions = self.setPositionList(self.codons)
        return positions[codon.upper()]


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
            Input: a set of amino acids in the DNA sequence
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
    
    def transcript(self,sequence):
        '''
            Input: a DNA sequence
            Output: a list of mRNAs
            Approach: Using Open Reading Frame (ORF) analysis
        '''
        # First, break down the sequence into codons
        codons = self.decipher(sequence.upper()) 
        
        # Second, define coding segments in the DNA sequence
        start_pos = self.getPositionList('ATG')
        stop_pos = []
        for pos,codon in codons.items():
            if codon in self.stop:
                stop_pos.append(pos)
        # if no stop codon in the DNA sequence
        if stop_pos == []:
            stop_pos = [len(codons.values())]
        
        coding_regions = dict(zip(start_pos,stop_pos))
        
        # Third, retrieve coding DNA sequence from start and stop positions
        pre_mRNAs = []
        for start,stop in coding_regions.items():
            pre_mRNA = self.getSubsequence(start,stop)
            pre_mRNA = pre_mRNA.replace("T","U")
            pre_mRNAs.append(pre_mRNA)
        
        mRNAs = "".join(pre_mRNAs) # this is biologically incorrect because we need to remove non-coding regions within each pre_mRNA before generating mRNA
        
        return mRNAs

    def translate(self,mRNA):
        codons = self.decipher(mRNA)
        amino_acids = self.getAminoAcidName('rna',codons.values())
        for codon,aa in amino_acids.items():
            if aa == 'stop':
                amino_acids[codon] = "*"
            elif aa == 'Incomplete':
                amino_acids[codon] = ''
            else:
                amino_acids[codon]= aa[0].upper()
        
        aa = []
        for codon in codons.values():
            aa.append(amino_acids[codon])

        proteins ="".join(aa)
        return proteins

    # def labelMutation(self,mutation):
    #     mutation_type = ["synonymous","missense","nonsense"]

    