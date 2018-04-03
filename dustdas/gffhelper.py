#!/usr/bin/env python
import re
import sys
import json

import dustdas.fastahelper as fh


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
                   "attribute": gffcols[8].strip(),
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
        self.fasta_header = None
        self.fasta_sequence = None

    def to_json(self, omit_fasta=False):
        if omit_fasta:
            res = dict()
            for k,v in self.__dict__.items():
                if k in ["fasta_header", "fasta_sequence"]:
                    pass
                else:
                    res[k] = v
            return json.dumps(res, default=lambda o:o.__dict__, sort_keys=True, indent=4)

        return json.dumps(self, default=lambda o: o.__dict__,
                          sort_keys=True, indent=4)
        #return json.dumps(self) #todo


    def get_sequence(self, fastafile=None, fastadct=None, regex=None):
        if fastafile:
            for header, seq in fh.FastaParser.read_fasta(fasta=fastafile):
                if regex:
                    p = re.compile(regex)
                    m = p.match(header)
                    if m:
                        return header, seq
                else:
                    if header == self.seqname:
                        return header, seq
        elif fastadct:
            for header, seq in fastadct.items():
                if regex:
                    p = re.compile(regex)
                    m = p.match(header)
                    if m:
                        return header, seq
                else:
                    return header, fastadct[header]





    def attrib_filter_fun(self, tfun=None, targ=None, vfun=lambda x,y : x.startswith(y), varg=None):
        if tfun and vfun:
            for a in self.attributes:
                t = tfun(a.tag, targ)
                v = vfun(a.value, varg)
                if t and v:
                    return True

        ##if tfun and targ:
        # #  return  tfun(self.tag, targ)
        #if vfun and varg:
        #    return vfun(self.value, varg)
        return False #TODO


    def attrib_filter(self, tag=None, value=None): # todo add passing filter with arbitrary function
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
        self.fasta_header = header
        self.fasta_sequence = seq


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
