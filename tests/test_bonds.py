#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pyformulate.smiles import loads
from pyformulate.smiles.models import CisTransType


ELEMENTS = ("F", "Cl", "Br", "I")


def test_implicit_single_bonds():
    for c in ELEMENTS:
        assert loads(2 * f"[{c}]").molecules[0].atoms()[0].valency == 1


def test_explicit_single_bonds():
    for c in ELEMENTS:
        assert loads(f"[{c}]-[{c}]").molecules[0].atoms()[0].valency == 1


def test_double_bonds():
    for c in ELEMENTS:
        assert loads(f"[{c}]=[{c}]").molecules[0].atoms()[0].valency == 2


def test_triple_bonds():
    for c in ELEMENTS:
        assert loads(f"[{c}]#[{c}]").molecules[0].atoms()[0].valency == 3


def test_quadruple_bonds():
    for c in ELEMENTS:
        assert loads(f"[{c}]$[{c}]").molecules[0].atoms()[0].valency == 4


def test_cis_trans_bonds():
    for c in ELEMENTS:
        bond_types = list(loads(f"[{c}]/[{c}]").molecules[0].atoms()[0].bonds.values())
        assert bond_types[0].cis_trans == CisTransType.BOTTOM_TOP
