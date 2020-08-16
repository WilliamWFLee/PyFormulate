#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pyformulate.smiles import loads
from pyformulate.smiles.models import ChiralClass


def test_chiral_classes():
    for chiral_class in ChiralClass:
        if not chiral_class.name.startswith("_"):
            assert (
                loads("[C@{0.name}]".format(chiral_class))
                .molecules[0]
                .atoms()[0]
                .chiral_class
                == chiral_class
            )
