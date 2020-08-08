#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from enum import Enum
from typing import Tuple

from .atom import Atom


class BondType(Enum):
    """
    Bond type enumeration
    """

    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    QUADRUPLE = 4


class Bond:
    """
    Represents a bond between two atoms
    """

    def __init__(self, atoms: Tuple[Atom, Atom], type_: BondType):
        self.atoms = atoms
        self.type = type_
