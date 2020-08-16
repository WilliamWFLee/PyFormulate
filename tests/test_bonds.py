#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pyformulate.smiles import loads


def test_implicit_single_bonds():
    for c in ("F", "Cl", "Br", "I"):
        assert (
            sum(1 for bond in loads(2 * c).molecules[0].bonds if bond.order == 1) == 1
        )


def test_explicit_single_bonds():
    for c in ("F", "Cl", "Br", "I"):
        assert (
            sum(1 for bond in loads(f"{c}-{c}").molecules[0].bonds if bond.order == 1)
            == 1
        )


def test_double_bonds():
    for c in "CNO":
        assert (
            sum(1 for bond in loads(f"{c}={c}").molecules[0].bonds if bond.order == 2)
            == 1
        )


def test_triple_bonds():
    for c in "CN":
        assert (
            sum(1 for bond in loads(f"{c}={c}").molecules[0].bonds if bond.order == 2)
            == 1
        )


def test_quadruple_bonds():
    assert sum(1 for bond in loads("C$C").molecules[0].bonds if bond.order == 4) == 1
