import pytest
from dustdas import gffhelper
import os
dir = os.path.dirname(__file__)

@pytest.mark.parametrize("gff", [os.path.join(dir,'test.gff3')])
@pytest.mark.parametrize("expected", [['gff-version 3.2.1', 'sequence-region ctg123 1 1497228']])
def test_metadata(gff, expected):
    g = gffhelper.GFFFile(gff)
    assert g.path == gff
    assert g.metadata == expected


def test_seqid():
    pass


def test_source():
    pass


def test_type():
    pass


def test_start():
    pass


def test_end():
    pass


def test_score():
    pass


def test_strand():
    pass


def test_phase():
    pass


def test_attributes():
    pass

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

