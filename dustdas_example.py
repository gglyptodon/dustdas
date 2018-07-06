#!/usr/bin/env python
import re
import sys
import argparse
from dustdas import gffhelper as gh
from dustdas import fastahelper as fh


def format_help():
    print("""
        from http://www.ensembl.org/info/website/upload/gff.html:
    
    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature. 
    #----------------------------------------------------------#

    """)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", type=str, help="path to gff file (example works with phytozome data)")
    parser.add_argument("fasta", type=str, help="path to genome fasta file")
    args = parser.parse_args()
    allthefeatures = {}

    genes = []
    five_prime_utrs = []
    three_prime_utrs = []
    mrnas = []
    exons = []
    cds = []

    for o in gh.read_gff_file(infile=args.gff):
        if o.feature in allthefeatures:
            allthefeatures[o.feature] += 1
        else:
            allthefeatures[o.feature] = 1

        if o.feature == "mRNA":
            mrnas.append(o)

        if o.feature == "gene":
            genes.append(o)

        if o.feature == "CDS":
            cds.append(o)

        if o.feature == "exon":
            exons.append(o)

        if o.feature == "five_prime_UTR":
            five_prime_utrs.append(o)

        if o.feature == "three_prime_UTR":
            three_prime_utrs.append(o)

    print(allthefeatures)

    # show first few genes
    for g in genes[0:2]:
        print("gene:", g)
        print("mRNA:")
        # print all mRNAs whose "ID" tag value of mrna starts with "Name" tag value of gene
        print([m for m in mrnas if
               m.attrib_filter_fun(tfun=lambda x, y: x == y, targ="ID", vfun=lambda x, y: x.startswith(y),
                                   varg=[a.value for a in g.attributes if a.tag == "Name"][0])])
        print("exons:")
        print([(e.start, e.end, e.strand, e) for e in exons if
               e.attrib_filter_fun(tfun=lambda x, y: x == y, targ="ID", vfun=lambda x, y: x.startswith(y),
                                   varg=[a.value for a in g.attributes if a.tag == "Name"][0])])
        print("cds:")
        print([(e.start, e.end, e.strand) for e in cds if
               e.attrib_filter_fun(tfun=lambda x, y: x == y, targ="ID", vfun=lambda x, y: x.startswith(y),
                                   varg=[a.value for a in g.attributes if a.tag == "Name"][0])])
        print("3pUTR:")
        print([(e.start, e.end, e.strand) for e in three_prime_utrs if
               e.attrib_filter_fun(tfun=lambda x, y: x == y, targ="ID", vfun=lambda x, y: x.startswith(y),
                                   varg=[a.value for a in g.attributes if a.tag == "Name"][0])])
        print("5pUTR:")
        print([(e.start, e.end, e.strand) for e in five_prime_utrs if
               e.attrib_filter_fun(tfun=lambda x, y: x == y, targ="ID", vfun=lambda x, y: x.startswith(y),
                                   varg=[a.value for a in g.attributes if a.tag == "Name"][0])])

        print("#########")

    fasta_dict = fh.FastaParser.read_fasta_whole(args.fasta)

    def setupseq(gffobj, fastadct, regex):

        """

        :param fastadct: dict
        :param regex: rawstring for regex
        :type gffobj: GFFObject
        """
        ident = gffobj.seqname
        r = regex.format(ident)
        h, s = gffobj.get_sequence(fastadct=fastadct, regex=r)
        seq = fh.FastaParser.get_sequence_by_coordinates(orig_seq=s,
                                                         start=gffobj.start,
                                                         end=gffobj.end,
                                                         strand=gffobj.strand)
        try:
            gffobj.attach_fasta(";".join([h, gffobj.start,
                                          gffobj.end,
                                          gffobj.strand,
                                          gffobj.frame]),
                                seq)
        except fh.SequenceTranslationException as e:
            print("{}, skipping {}".format(e, gffobj), file=sys.stderr)

    for c in cds:
        setupseq(c, fasta_dict, r"^{}")
    for m in mrnas[:11]:
        identifier = [a.value for a in m.attributes if a.tag == "ID"][0]
        # find all cds that belong to mrna
        conc = ""
        with open("{}_cds.json".format(identifier), 'w') as out:
            mrna_cds = [c for c in cds if c.attrib_filter_fun(tfun=lambda x, y: x == y,
                                                              targ="Parent",
                                                              vfun=lambda x, y: x.startswith(y),
                                                              varg=[a.value for a in m.attributes if a.tag == "ID"][0])
                        ]
            for c in mrna_cds:
                out.write(c.to_json())
                conc += c.fasta_sequence
        with open("{}_cds_concatenated.fasta".format(identifier), 'w') as outfa:
            outfa.write(">{}\n{}\n".format(identifier + "_cds", conc))
        with open("{}_cds_prot_concatenated.fasta".format(identifier), 'w') as outfa:
            fshift = 0
            if m.frame != ".":
                fshift = m.frame
            outfa.write("{}\n{}".format(identifier + "_cds_protein", fh.SeqTranslator.dna2prot(conc, frameshift=0)))


if __name__ == "__main__":
    main()
