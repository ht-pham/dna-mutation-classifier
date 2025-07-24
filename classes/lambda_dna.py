from classes.dna import DNA

class LambdaGenome(DNA):

    def __init__(self,name,dna_sequence):
        super().__init__()
        self.name = name
        self.sequence = dna_sequence
        self.amino_acids = set()
    def __str__(self):
        return super().__str__()

    def load_from_file(self):
        from Bio import SeqIO
        record = SeqIO.parse("data/lambda_genome.fasta", "fasta")
        self.sequence = next(record).seq
        
    def load_genome(self,file_path="data/lambda_genome.fasta"):
        return super().load_genome(file_path)

    def getSequence(self):
        _,self.sequence = super().load_genome("data/lambda_genome.fasta")
        return self.sequence
    
    def generate(self, num_nucleotides):
        return super().generate(num_nucleotides)
    
    def decipher(self,sequence=getSequence):
        return super().decipher(sequence)
       

    def setPositionList(self, codons_dict):
        return super().setPositionList(codons_dict)
    
    def getPositionList(self, codon):
        return super().getPositionList(codon)

    def getAminoAcidAt(self, codon_index):
        return super().getAminoAcidAt(codon_index)
    
    def getAminoAcidName(self,nucleoic_acid,amino_acids):
        return super().getAminoAcidName(nucleoic_acid,amino_acids)
    
    def countAll(self,nucleoic_acid,sequence):
        return super().countAll(nucleoic_acid,sequence)
    
    def printStatistics(self,nucleoic_acid='dna',sequence=getSequence):
        # codons = self.decipher(sequence)
        # #print(f"Description: {self.name}")
        # print(f"Number of codons: {len(codons)}")
        # print(f"Number of amino acids: {len(set(codons.values()))}")
        return super().printStatistics(nucleoic_acid,sequence)
    
    def transcribe(self,sequence):
        return super().transcribe(sequence)
    
    def alterSNP(self, original_sequence='ATGTAG'):
        return super().alterSNP(original_sequence)