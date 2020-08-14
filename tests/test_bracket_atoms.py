#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest

from pyformulate.smiles import loads
from pyformulate.smiles.decoder import DecodeError, DecodeWarning
from pyformulate.smiles.models import Element


def test_known_elements():
    for elem in Element:
        if elem == Element.UNKNOWN:
            continue
        assert len(loads(f"[{elem.name}]").molecules[0].atoms(elem)) == 1


def test_unknown_element():
    assert len(loads("[*]").molecules[0].atoms(Element.UNKNOWN)) == 1


def test_invalid_elements():
    with pytest.raises(DecodeError):
        loads("[]")  # Empty bracket atom

    # Invalid elements
    for element in ("Pq", "Ny", "Zh", "Pi", "Nj", "Af"):
        with pytest.raises(DecodeError):
            loads(f"[{element}]")

    # Invalid aromatics
    for element in ("sd", "t", "v", "sd", "ap", "sg"):
        with pytest.raises(DecodeError):
            loads(f"[{element}]")


def test_explicit_hydrogen_count():
    assert loads("[CH4]").molecules[0].elem_count(Element.H) == 4
    assert loads("[ClH1]").molecules[0].elem_count(Element.H) == 1


def test_implicit_hydrogen_count():
    assert loads("[ClH]").molecules[0].elem_count(Element.H) == 1


def test_explicit_no_hydrogen_count():
    assert loads("[ClH0]").molecules[0].elem_count(Element.H) == 0


def test_implicit_no_hydrogen_count():
    assert loads("[Cl]").molecules[0].elem_count(Element.H) == 0


def test_illegal_hydrogen_count():
    with pytest.raises(DecodeError):
        loads("[HH1]")


def test_no_charge():
    assert loads("[Cl]").molecules[0].atoms()[0].charge == 0


def test_single_charge():
    assert loads("[H+]").molecules[0].atoms()[0].charge == 1
    assert loads("[Cl-]").molecules[0].atoms()[0].charge == -1


def test_double_charge():
    with pytest.warns(DecodeWarning):
        assert loads("[Ca++]").molecules[0].atoms()[0].charge == 2
        assert loads("[O--]").molecules[0].atoms()[0].charge == -2


def test_multiple_charge():
    assert loads("[Cu+3]").molecules[0].atoms()[0].charge == 3
    assert loads("[Ca+2]").molecules[0].atoms()[0].charge == 2
    assert loads("[Cr+6]").molecules[0].atoms()[0].charge == 6
    assert loads("[S-2]").molecules[0].atoms()[0].charge == -2


def test_isotope():
    assert loads("[13CH4]").molecules[0].atoms(Element.C)[0].isotope == 13
    assert loads("[2H+]").molecules[0].atoms()[0].isotope == 2
    assert loads("[238U]").molecules[0].atoms()[0].isotope == 238


def test_trailing_zeroes_isotope():
    assert loads("[002CH4]").molecules[0].atoms(Element.C)[0].isotope == 2


def test_zero_isotope():
    assert loads("[0Ca]").molecules[0].atoms()[0].isotope == 0
