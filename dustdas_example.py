#!/usr/bin/env python
import re
import sys
import argparse
from dustdas import gffhelper as gh
from dustdas import fastahelper as fh


def formatHelp():
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
    five_prime_UTRs = []
    three_prime_UTRs = []
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
            five_prime_UTRs.append(o)

        if o.feature == "three_prime_UTR":
            three_prime_UTRs.append(o)


    print(allthefeatures)

    # show first few genes
    for g in genes[0:2]:
        print("gene:", g)
        print("mRNA:")
        # print all mRNAs whose "ID" tag value of mrna starts with "Name" tag value of gene
        print ([m for m in mrnas if m.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in g.attributes if a.tag=="Name"][0])])
        print("exons:")
        print ([(e.start,e.end,e.strand,e)  for e in exons if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in g.attributes if a.tag=="Name"][0])])
        print("cds:")
        print ([(e.start,e.end,e.strand)  for e in cds if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in g.attributes if a.tag=="Name"][0])])
        print("3pUTR:")
        print ([(e.start,e.end,e.strand)  for e in three_prime_UTRs if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in g.attributes if a.tag=="Name"][0])])
        print("5pUTR:")
        print ([(e.start,e.end,e.strand)  for e in five_prime_UTRs if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in g.attributes if a.tag=="Name"][0])])

        print ("#########")

    fasta_dict = fh.FastaParser.read_fasta_whole(args.fasta)
    for m in mrnas[0:3]:
        identifier = m.seqname # this will be used to match the (genome) fasta sequence
        # for cds or protein fasta this could be sth like: identifier = m.attrib_filter(tag="pacid")
        # and regex pattern would be sth like this: r = r"""(.*)pacid={}""".format(identifier[0].value)

        # matching fasta header starts with identifier
        r = r"""^{}""".format(identifier)

        # attach fasta sequence to object
        h,s = m.get_sequence(fastadct=fasta_dict, regex=r)
        mrnaseq=fh.LargeFastaParser().getSequenceByCoordinatesOnLandmark(landmarkSeq=s,start=m.start, end=m.end, strand=m.strand)
        m.attach_fasta(";".join([h,m.start, m.end, m.strand]),mrnaseq)


        # look at all exons of mrna
        # filter for all exons whose parent tag value starts with mrna ID tag
        mrnaExons = [e for e in exons if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "Parent", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in m.attributes if a.tag=="ID"][0])]

        for e in mrnaExons:
            print(e)
            #print(fh.LargeFastaParser().getSequenceByCoordinatesOnLandmark(landmarkSeq=m.fasta_sequence, start=e.start,
            #                                                          end=e.end, strand=e.strand))

        print ("#")

    with open ("tmp.json", 'w') as out:
        for m in exons[0:3]:
            out.write(m.to_json())
        for m in mrnas[0:3]:
            out.write(m.to_json(omit_fasta=True))
        for m in mrnas[0:3]:
            out.write(m.to_json(omit_fasta=False))







if __name__ == "__main__":
    main()
