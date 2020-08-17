#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import string

import pytest

from pyformulate.models import Element
from pyformulate.smiles import loads
from pyformulate.smiles.decoder import ALIPHATIC_ORGANIC, DecodeError


def test_valid_aliphatics():
    for elem in ALIPHATIC_ORGANIC:
        loads(elem)


def test_invalid_aliphatics():
    for elem in "".join(
        filter(lambda x: x not in ALIPHATIC_ORGANIC, string.ascii_uppercase)
    ):
        with pytest.raises(DecodeError):
            loads(elem)


def test_valid_aromatics():
    loads("c1ccccc1")  # Benzene
    loads("c1cocc1")  # Furan


def test_implicit_hydrogens():
    assert loads("C").molecules[0].elem_count(Element.H) == 4
    assert loads("N").molecules[0].elem_count(Element.H) == 3
    assert loads("O").molecules[0].elem_count(Element.H) == 2
