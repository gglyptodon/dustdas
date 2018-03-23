#!/usr/bin/env python
import re
import sys
import getopt
import argparse


class GFFObject(object):
    @staticmethod
    def parse_gffline(gffline):
        """takes one line of a gff3 file and returns a dictionary"""

        if gffline.startswith("#"):
            pass
        else:
            gffcols = gffline.split("\t")

            res = {"seqname": gffcols[0],
                   "source": gffcols[1],
                   "feature": gffcols[2],
                   "start": gffcols[3],
                   "end": gffcols[4],
                   "score": gffcols[5],
                   "strand": gffcols[6],
                   "frame": gffcols[7],
                   "attribute": gffcols[8],
                   }

            return res

    def __init__(self, gffline):
        d = GFFObject.parse_gffline(gffline)
        self.seqname = d["seqname"]
        self.source = d["source"]
        self.feature = d["feature"]
        self.start = d["start"]
        self.end = d["end"]
        self.score = d["score"]
        self.strand = d["strand"]
        self.frame = d["frame"]
        self.attribute = [GFFAttribute(x.strip()) for x in d["attribute"].split(";")]



class GFFAttribute(object):
    def __init__(self, attribute_str):
        #attribute_list = attribute_str.split(";")
        p = re.compile(r"""(.*)=(.*)""")
        res = []
        #for a in attribute_list:
        m = p.match(attribute_str)
        if m:
            r = [{"tag": m.groups()[0], "value": m.groups()[1]}]
            self.tag = m.groups()[0]
            self.value = m.groups()[1]
        else:
            self.tag = "wat"
            self.value = "wat"

    def __repr__(self):
        return "tag: {} - value: {}".format(self.tag, self.value)

def usage():
    """Print usage."""

    sys.exit(2)


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


def checkGff3(infile):
    with open(infile, 'r') as f:
        l = f.readline().strip()
        if l == "##gff-version 3":
            return True
        else:
            return False


def read_gff_attrib(infile, attrib):
    with open(infile, 'r') as f:
        for l in f:
            if l.strip() == "":
                pass
            elif l.startswith("#"):
                pass
            else:
                tmp = l.split("\t")[2].lower().strip()
                if tmp == attrib.lower():
                    ln = l.split("\t")[8].strip()
                    yield (ln)

def filtersth(tag, value):
    pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", type=str, help="path to gff file")
    args = parser.parse_args()
    with open(args.gff) as inf:
        for l in inf:
            if l.startswith("#"):
                pass
            else:
                g = GFFObject(gffline=l)
                #print(g.attribute)
                # show all available attributes
                if "Parent" in [x.tag for x in g.attribute]:
                    print ([a.value for a in g.attribute if a.tag == "Parent"])
                #

if __name__ == "__main__":
    main()
