#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pyformulate.smiles import loads


def test_separate_molecules():
    assert len(loads("[Na+].[Cl-]").molecules) == 2
    assert len(loads("C.C.C").molecules) == 3


def test_connected_molecules():
    assert len(loads("C1.C1").molecules) == 1
    assert len(loads("C12.C1.C2").molecules) == 1
    assert len(loads("C1(.C1)").molecules) == 1
    assert len(loads("C12(.C1(.C2))").molecules) == 1


def test_dot_bond_branches():
    assert len(loads("C(.C)").molecules) == 2
    assert len(loads("C(.C(.C))").molecules) == 3
