#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Models specific to SMILES interpretation
"""

from enum import Enum
from typing import Optional

from .. import models
from ..models import CisTransType, Element

# Defines the organic subset of elements
ALIPHATIC_ORGANIC = ("B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I")
AROMATIC_ORGANIC = "bcnosp"

# Defines the aromatics
AROMATIC = ("b", "c", "n", "o", "p", "s", "se", "as")

CHAR_TO_BOND_ORDER = {
    "-": 1,
    "=": 2,
    "#": 3,
    "$": 4,
    ":": 1,
    "/": 1,
    "\\": 1,
}


CHAR_TO_CIS_TRANS = {
    "\\": CisTransType.TOP_BOTTOM,
    "/": CisTransType.BOTTOM_TOP,
}

# The "normal" valencies of atoms, as defined in SMILES
VALENCIES = {
    Element.B: (3,),
    Element.C: (4,),
    Element.N: (3, 5),
    Element.O: (2,),
    Element.P: (3, 5),
    Element.S: (2, 4, 6),
    Element.F: (1,),
    Element.Cl: (1,),
    Element.Br: (1,),
    Element.I: (1,),
    Element.At: (1,),
    Element.Ts: (1,),
}


class ChiralClass(Enum):
    """
    Enumeration for the tags used in SMILES for representing chiral stereochemistry
    """

    _ANTICLOCKWISE = -1  # Indicates an atom requires further analysis
    _CLOCKWISE = -2  # to determine chiral class
    TH1 = 1
    TH2 = 2
    AL1 = 3
    AL2 = 4
    SP1 = 5
    SP2 = 6
    SP3 = 7
    TB1 = 8
    TB2 = 9
    TB3 = 10
    TB4 = 11
    TB5 = 12
    TB6 = 13
    TB7 = 14
    TB8 = 15
    TB9 = 16
    TB10 = 17
    TB11 = 18
    TB12 = 19
    TB13 = 20
    TB14 = 21
    TB15 = 22
    TB16 = 23
    TB17 = 24
    TB18 = 25
    TB19 = 26
    TB20 = 27
    OH1 = 28
    OH2 = 29
    OH3 = 30
    OH4 = 31
    OH5 = 32
    OH6 = 33
    OH7 = 34
    OH8 = 35
    OH9 = 36
    OH10 = 37
    OH11 = 38
    OH12 = 39
    OH13 = 40
    OH14 = 41
    OH15 = 42
    OH16 = 43
    OH17 = 44
    OH18 = 45
    OH19 = 46
    OH20 = 47
    OH21 = 48
    OH22 = 49
    OH23 = 50
    OH24 = 51
    OH25 = 52
    OH26 = 53
    OH27 = 54
    OH28 = 55
    OH29 = 56
    OH30 = 57


class Atom(models.Atom):
    def __init__(
        self,
        element: Element,
        *,
        isotope: Optional[int] = None,
        charge: int = 0,
        aromatic: bool = False,
        atom_class: int = 0,
        chiral_class: Optional[ChiralClass] = None,
    ):
        super().__init__(
            element,
            isotope=isotope,
            charge=charge,
            aromatic=aromatic,
            atom_class=atom_class,
            chiral_class=chiral_class,
        )

    def __repr__(self):
        return (
            "{0.__qualname__}(element={1.element}, isotope={1.isotope}, "
            "charge={1.charge}, chiral_class={1.chiral_class}, "
            "aromatic={1.aromatic}, atom_class={1.atom_class})"
        ).format(type(self), self)


class Molecule(models.Molecule):
    pass
