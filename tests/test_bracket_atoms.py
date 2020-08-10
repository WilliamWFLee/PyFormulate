#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pyformulate.smiles import loads
from pyformulate.smiles.models import Element


def test_known_elements():
    for elem in Element:
        if elem == Element.UNKNOWN:
            continue
        assert loads(f"[{elem.name}]").molecules[0].atoms[0].element == elem


def test_unknown_element():
    assert loads("[*]").molecules[0].atoms[0].element == Element.UNKNOWN
