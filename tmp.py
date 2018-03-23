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
    for o in gh.read_gff_file(infile=args.gff):
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
    print (allthefeatures)

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
