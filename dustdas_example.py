#!/usr/bin/env python
import re
import sys
import argparse
import json

from dustdas import gffhelper as gh
from dustdas import fastahelper as fh


def format_help():
    print("""
    1 - seqid:  ID of the landmark used to establish the coordinate system for the current feature.
    2 - source: data source or generating software
    3 - type:   a term or accession from the SOFA sequence ontology
    4 - start:  int, start position, numbering starting at 1
    5 - end:    int, end position, numbering starting at 1
    6 - score:  floating point score
    7 - strand: '+' for positive strand (relative to the landmark),
            '-' for minus strand, and . for features that are not stranded, 
            '?' for relevant but unknown
    8 - phase: '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon etc
    9 - attributes: semicolon-separated list of tag-value pairs, e.g.  
        ID, Name, Alias, Parent, Target, Gap, Derives_from, Note, Dbxref, Ontology_term, Is_circular, ...
        "All attributes that begin with an uppercase letter are reserved for later use. 
        Attributes that begin with a lowercase letter can be used freely by applications." 
        (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
        
    See also:
        https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        http://www.ensembl.org/info/website/upload/gff3.html
    #----------------------------------------------------------#

    """)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", type=str, help="path to gff file (example works with phytozome data)")
    parser.add_argument("fasta", type=str, help="path to genome fasta file")
    args = parser.parse_args()
    all_features = {}

    genes = []
    five_prime_utrs = []
    three_prime_utrs = []
    mrnas = []
    exons = []
    cds = []

    for o in gh.read_gff_file(infile=args.gff):
        if o.type in all_features:
            all_features[o.type] += 1
        else:
            all_features[o.type] = 1

        if o.type == "mRNA":
            mrnas.append(o)

        if o.type == "gene":
            genes.append(o)

        if o.type == "CDS":
            cds.append(o)

        if o.type == "exon":
            exons.append(o)

        if o.type == "five_prime_UTR":
            five_prime_utrs.append(o)

        if o.type == "three_prime_UTR":
            three_prime_utrs.append(o)

    print(all_features)

    # show first few genes
    #for g in genes[0:2]:
    #    print("gene:", g)
    #    print("mRNA:")
    #    print (g)

    fasta_dict = fh.FastaParser.read_fasta_whole(args.fasta)

    def setupseq(gffobj, fastadct, regex):
        """
        :param fastadct: dict
        :param regex: rawstring for regex
        :type gffobj: GFFObject
        """
        ident = gffobj.seqid
        r = regex.format(ident) # eg "^{} .*"
        h, s = gffobj.get_sequence(fastadct=fastadct, regex=r)
        seq = fh.FastaParser.get_sequence_by_coordinates(orig_seq=s,
                                                         start=gffobj.start,
                                                         end=gffobj.end,
                                                         strand=gffobj.strand)
        try:
            gffobj.attach_fasta(";".join([h, gffobj.start,
                                          gffobj.end,
                                          gffobj.strand,
                                          gffobj.phase]),
                                seq)
        except fh.SequenceTranslationException as e:
            print("{}, skipping {}".format(e, gffobj), file=sys.stderr)


    #example: exons for first mrnas
    for m in mrnas[4:6]:
            with open ("{}_exons.json".format(m.get_ID()),'w') as out:
                #out.write("[")
                for e in [x for x in exons if  m.get_ID() in x.get_Parent()]:
                    setupseq(e, fasta_dict, r"^{} .*")
                    #out.write(e.to_json())
                out.write(json.dumps([x for x in exons if  m.get_ID() in x.get_Parent()], default=lambda o: o.__dict__, sort_keys=True, indent=4))
               # out.write("]")
if __name__ == "__main__":
    main()
