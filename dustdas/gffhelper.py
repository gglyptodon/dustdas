#!/usr/bin/env python
import re
import sys

import  dustdas.fastahelper  as fh

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
        self.attribute = d["attribute"]
        self.attributes = [GFFAttribute(x.strip()) for x in d["attribute"].split(";")]

    def get_sequence(self, fastafile, regex=None):
        """ example:(.*)ID=AHYPO_002876-RA.v1.0 """
        #header, seq = fh.FastaParser().read_fasta(fasta=fastafile)
        for header, seq in fh.FastaParser().read_fasta(fasta=fastafile):
            #print(header)
            if regex:
                p = re.compile(regex)
                m = p.match(header)
                if m:
                    return header, seq
            else:
                if header == self.seqname:
                    return header, seq

    def attrib_filter(self, tag=None, value=None):
        if tag and not value:
            return ([a for a in self.attributes if a.tag == tag])
        elif value and not tag:
            return ([a for a in self.attributes if a.value == value])
        elif value and tag:
            return ([a for a in self.attributes if a.value == value and a.tag == tag])
        else:
            print("needs tag or value to filter. returns list of matches", file=sys.stderr)

    def __repr__(self):
        return "{},{},{},{},{},{},{},{},{}".format(self.seqname, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, self.attributes)

    def attach_fasta(self, header, seq):
        self.fastaheader = header
        self.fastaseqence = seq


class GFFAttribute(object):
    def __init__(self, attribute_str):
        p = re.compile(r"""(.*)=(.*)""")
        m = p.match(attribute_str)
        if m:
            r = [{"tag": m.groups()[0], "value": m.groups()[1]}]
            self.tag = m.groups()[0]
            self.value = m.groups()[1]
        else:
            self.tag = "wat"
            self.value = "wat"

    def __repr__(self):
        return "<tag:{},value:{}>".format(self.tag, self.value)


def read_gff_file(infile):
    with open(infile, 'r') as f:
        for l in f:
            if l.strip() == "":
                pass
            elif l.startswith("#"):
                pass
            else:
                obj = GFFObject(gffline=l)
                yield obj
                #tmp = l.split("\t")[2].lower().strip()
                #if tmp == attrib.lower():
                #    ln = l.split("\t")[8].strip()
                #    yield (ln)

def make_filter(tag, value):
    pass
