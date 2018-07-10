import pytest
from dustdas import gffhelper
import os
dir = os.path.dirname(__file__)


@pytest.mark.parametrize("gff, expected", [
    (os.path.join(dir,'test.gff3'),['gff-version 3.2.1', 'sequence-region ctg123 1 1497228']),
    (os.path.join(dir,'test2.gff3'),['gff-version 3.2.1', 'sequence-region ctg1234 1 1497228'])
    ])
def test_metadata(gff, expected):
    g = gffhelper.GFFFile(gff)
    assert g.path == gff
    assert g.metadata == expected


@pytest.mark.parametrize("gff, expected", [
    (os.path.join(dir,'test.gff3'),'ctg123'),
    (os.path.join(dir,'test2.gff3'),'ctg1234')
    ])
def test_seqid(gff, expected):
    g = gffhelper.GFFFile(gff)
    for o in g.get_gff_objects():
        assert o.seqid == expected

@pytest.mark.parametrize("gff, expected", [
    (os.path.join(dir,'test.gff3'),'.'),
    (os.path.join(dir,'test2.gff3'),'.')
    ])
def test_source(gff, expected):
    g = gffhelper.GFFFile(gff)
    for o in g.get_gff_objects():
        assert o.source == expected


@pytest.mark.parametrize("gff, index, expected", [
    (os.path.join(dir,'test.gff3'),0,'gene'),
    (os.path.join(dir,'test.gff3'),1, 'TF_binding_site'),
    (os.path.join(dir,'test.gff3'),2, 'mRNA'),
    (os.path.join(dir,'test.gff3'),4, 'mRNA'),
    (os.path.join(dir,'test.gff3'),22, 'CDS'),
    (os.path.join(dir,'test3.gff3'),0, 'mRNA'),
    (os.path.join(dir,'test3.gff3'),1, 'exon')
    ])
def test_type(gff, index, expected):
    g = gffhelper.GFFFile(gff)
    objs = list(g.get_gff_objects())
    assert objs[index].type == expected

@pytest.mark.parametrize("gff, index, expected_start, expected_end", [
    (os.path.join(dir,'test.gff3'),0, 1000, 9000),
    (os.path.join(dir,'test.gff3'),1, 1000, 1012),
    (os.path.join(dir,'test.gff3'),2, 1050, 9000),
    (os.path.join(dir,'test.gff3'),4, 1300, 9000),
    (os.path.join(dir,'test.gff3'),22, 7000, 7600),
    (os.path.join(dir,'test3.gff3'),0, 1300, 9000),
    (os.path.join(dir,'test3.gff3'),1, 1300, 1500)
    ])
def test_start_end(gff, index, expected_start, expected_end):
    g = gffhelper.GFFFile(gff)
    objs = list(g.get_gff_objects())
    assert int(objs[index].start) == expected_start
    assert int(objs[index].end) == expected_end


@pytest.mark.parametrize("gff, index, expected_score, expected_strand, expected_phase", [
    (os.path.join(dir,'test.gff3'),0, ".","+", "."),
    (os.path.join(dir,'test.gff3'),1, ".","+", "."),
    (os.path.join(dir,'test.gff3'),2, ".","+", "."),
    (os.path.join(dir,'test.gff3'),4, ".","+", "."),
    (os.path.join(dir,'test.gff3'),20, ".","+", "0"),
    (os.path.join(dir,'test.gff3'),22, ".","+", "1"),
    (os.path.join(dir,'test3.gff3'),1, 1e-2,"+", "."),
    (os.path.join(dir,'test3.gff3'),2, 0.9, "+", "."),
    (os.path.join(dir,'test3.gff3'),3, 0.000001, "+", ".")
    ])
def test_score(gff, index, expected_score, expected_strand, expected_phase):
    g = gffhelper.GFFFile(gff)
    objs = list(g.get_gff_objects())
    assert objs[index].score == expected_score
    assert objs[index].strand == expected_strand
    assert objs[index].phase == expected_phase


@pytest.mark.parametrize("gff, index, expected_attributes", [
    (os.path.join(dir,'test.gff3'),0, [{"tag":"ID", "value":"gene00001"},{"tag":"Name" , "value":"EDEN"}]),
    (os.path.join(dir,'test.gff3'),1, [{"tag":"ID", "value":"tfbs00001"},{"tag":"Parent" , "value":"gene00001"}]),
    (os.path.join(dir,'test.gff3'),6, [{"tag":"ID", "value":"exon00002"},{"tag":"Parent" , "value":"mRNA00001,mRNA00002"}]),
    (os.path.join(dir,'test.gff3'),9, [{"tag":"ID", "value":"exon00005"},{"tag":"Parent" , "value":"mRNA00001,mRNA00002,mRNA00003"}])
])
def test_attributes(gff, index, expected_attributes ):
    g = gffhelper.GFFFile(gff)
    objs = list(g.get_gff_objects())
    for e in expected_attributes:
        vals = [a.value for a in objs[index].attributes if a.tag == e["tag"]]
        assert vals == [list(e["value"].split(","))] # todo maybe get rid of [[]]


# attributes
# parent, alias, note, dbxref, ontology_term attributes can have multiple values separated by comma
def test_attrib_id():
    pass


def test_attrib_name():
    pass


def test_attrib_alias():
    pass
def test_attrib_alias_multi():
    pass

def test_attrib_parent():
    pass


def test_attrib_parent_multi():
    pass


def test_attrib_target():
    pass


def test_attrib_gap():
    pass


def test_attrib_derives_from():
    pass


def test_attrib_note():
    pass

def test_attrib_note_multi():
    pass


def test_attrib_ontology_term():
    pass

def test_attrib_ontology_term_multi():
    pass

def test_attrib_dbxref():
    pass


def test_attrib_dbxref_multi():
    pass


def test_attrib_is_circular():
    pass





def test_mrna_cds():
    """find all cds that belong to mrna"""
    pass

def test_features_phytozome():
    """ find all features in gff file (as provided by phytozome)"""
    pass

def test_features_ensembl():
    """ find all features in gff file (as provided by ensembl)"""
    pass

def test_proteinconversion():
    """"""
    pass

