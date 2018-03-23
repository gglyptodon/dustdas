#!/usr/bin/env python
import re
import sys
import argparse
from dustdas import gffhelper as gh


def formatHelp():
    print("""
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
    parser.add_argument("gff", type=str, help="path to gff file")
    parser.add_argument("fasta", type=str, help="path to fasta file")
    args = parser.parse_args()
    allthefeatures = {}
    genes = []
    five_prime_UTRs = []
    three_prime_UTRs = []
    mrnas = []
    exons = []
    cds = []

    for o in gh.read_gff_file(infile=args.gff):


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


        #print(o)
        t = o.attrib_filter(tag="Parent", value="AHYPO_002876-RA.v1.0")
        t = o.attrib_filter(tag="pacid", value="32823782")
        if o.feature in allthefeatures:
            allthefeatures[o.feature]+=1
        else:
            allthefeatures[o.feature]= 1

        #if o.feature == "mRNA":
        t = o.attrib_filter(tag="pacid")
        if t:
            #print(o.seqname)
            r = r"""(.*)ID={}""".format(t[0].value)
            r = r"""(.*)pacid={}""".format(t[0].value)
            #print(">"+str(o))
            #h,s = o.get_sequence(fastafile=args.fasta, regex=r)
            #o.attach_fasta(header=h, seq=s)
            #print(o.fastaseqence)
    #print (allthefeatures)

    #print(genes[0])
    for g in genes[0:0]:
        print("gene:", g)
        print("mRNA:")
        print ([e for e in mrnas if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in g.attributes if a.tag=="Name"][0])])
        print("exons:")
        print ([(e.start,e.end,e.strand,e)  for e in exons if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in g.attributes if a.tag=="Name"][0])])
        print("cds:")
        print ([(e.start,e.end,e.strand)  for e in cds if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in g.attributes if a.tag=="Name"][0])])
        print("3pUTR:")
        print ([(e.start,e.end,e.strand)  for e in three_prime_UTRs if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in g.attributes if a.tag=="Name"][0])])
        print("5pUTR:")
        print ([(e.start,e.end,e.strand)  for e in five_prime_UTRs if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in g.attributes if a.tag=="Name"][0])])

#TODO:relative to mrna: gaps between exons , filter by mrna instead of gene

        print ("#########")

    for m in mrnas[0:2]:
        print("mrna", m)
        print("exons:")
        print ([(e.start,e.end,e.strand)  for e in exons if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "Parent", vfun = lambda x,y: x.startswith(y), varg=[a.value for a in m.attributes if a.tag=="ID"][0])])
        print ("#########")
    #print(exons[0:20])

    #print ([e for e in exons if e.attrib_filter_fun(tfun= lambda x,y: x==y, targ = "ID", vfun = lambda x,y: x.startswith(y), varg="AHYPO_000001-RA.v1.0.exon.")])
    """with open(args.gff) as inf:
        for l in inf:
            if l.startswith("#"):
                pass
            else:
                g = gh.GFFObject(gffline=l)
                #print(g.attribute)
                # show all available attributes
                if "Parent" in [x.tag for x in g.attribute]:
                    print ([a.value for a in g.attributes if a.tag == "Parent"])
    
                #
"""
if __name__ == "__main__":
    main()
