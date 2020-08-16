#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pyformulate.smiles import loads
from pyformulate.smiles.models import CisTransType


ELEMENTS = ("F", "Cl", "Br", "I")


def test_implicit_single_bonds():
    for c in ELEMENTS:
        assert loads(2 * "[{}]".format(c)).molecules[0].bonds[0].order == 1


def test_explicit_single_bonds():
    for c in ELEMENTS:
        assert loads("[{0}]-[{0}]".format(c)).molecules[0].bonds[0].order == 1


def test_double_bonds():
    for c in ELEMENTS:
        assert loads("[{0}]=[{0}]".format(c)).molecules[0].bonds[0].order == 2


def test_triple_bonds():
    for c in ELEMENTS:
        assert loads("[{0}]#[{0}]".format(c)).molecules[0].bonds[0].order == 3


def test_quadruple_bonds():
    for c in ELEMENTS:
        assert loads("[{0}]$[{0}]".format(c)).molecules[0].bonds[0].order == 4


def test_cis_trans_bonds():
    for c in ELEMENTS:
        assert (
            loads("[{0}]/[{0}]".format(c)).molecules[0].bonds[0].cis_trans
            == CisTransType.BOTTOM_TOP
        )
