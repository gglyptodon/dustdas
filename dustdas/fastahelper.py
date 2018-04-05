import re


class FastaHelper(object):
    @staticmethod
    def remove_newlines(s):
        yield s.replace("\n", "")

    @staticmethod
    def insert_newlines(s, every=80):
        return '\n'.join(s[i:i + every] for i in range(0, len(s), every))

    @staticmethod
    def complement(s):
        tr = str.maketrans('AGTCagtc', 'TCAGtcag')
        res = s.translate(tr)
        return res

    @staticmethod
    def reverse_complement(s):
        tr = str.maketrans('AGTCagtc', 'TCAGtcag')
        res = s.translate(tr)
        return res[::-1]


class FastaParser(object):
    @staticmethod
    def read_fasta(fasta):
        """
        read from fasta fasta file 'fasta'
        yield header, sequence
        """
        name = ""
        fasta = open(fasta, "r")
        while True:
            line = name or fasta.readline()
            if not line:
                break
            seq = []
            while True:
                name = fasta.readline()
                name = name.rstrip()
                if not name or name.startswith(">"):
                    break
                else:
                    seq.append(name)
            joinedSeq = "".join(seq)
            line = line[1:]
            yield (line.rstrip(), joinedSeq.rstrip())
        fasta.close()

    @staticmethod
    def read_fasta_whole(fasta):
        """
        read from fasta fasta file 'fasta'
        returns dict header:sequence
        """
        res = dict()
        for h, s in FastaParser.read_fasta(fasta=fasta):
            res[h] = s
        return res

    @staticmethod
    def get_sequence_by_coordinates(orig_seq, start, end, strand, no_reverse_complement=False):
        """
        return subsection from orig_seq between start and end (gff coordinates, starting at 1)
        if strand is -, return the reverse complement of the subsequence
        unless no_reverse_complement is set to True
        """
        start = int(start) - 1  # translate gff coordinates
        end = int(end)
        res = orig_seq[start:end]

        if strand == "-" and not no_reverse_complement:
            res = orig_seq[start:end]
            return FastaHelper.reverse_complement(res)
        if strand not in ["-","+","."]:
            raise Exception("strand in unknown orientation.")

        return res


class SeqTranslator(object):
    RNAmap = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
        "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
    DNAmap = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
    @staticmethod
    def triplets(s, frameshift=0):
        # chop string into 3 char slices
        for i in range(0, int(len(s) / 3)):
            yield s[i * 3 + frameshift:i * 3 + 3 + frameshift]

    @staticmethod
    def triplets1(s, frameshift=0):  # todo
        if frameshift == 3:
            frameshift = 0
        # chop string into 3 char slices
        for i in range(0, int(len(s) / 3)):
            if frameshift == 0:
                yield s[i * 3 + frameshift:i * 3 + 3 + frameshift]
            elif frameshift == 1:
                yield s[1 + i * 3 + frameshift:1 + i * 3 + 3 + frameshift]
            elif frameshift == 2:
                yield s[2 + i * 3 + frameshift:2 + i * 3 + 3 + frameshift]
            else:
                raise Exception("frameshift out of range")

    @staticmethod
    def dna2prot(s, frameshift=0):
        res = ""
        for a in SeqTranslator.triplets1(s, frameshift):
            try:
                res += SeqTranslator.DNAmap[a]
            except KeyError as e:
                raise Exception("{} not in DNAmap".format(a))
                print(e)
                res += ""
        return res

    @staticmethod
    def rna2prot(s, frameshift=0):
        res = ""
        for a in SeqTranslator.triplets1(s, frameshift):
            try:
                res += SeqTranslator.RNAmap[a]
            except KeyError as e:
                print(e)
                raise Exception("{} not in RNAmap".format(a))
        return res
