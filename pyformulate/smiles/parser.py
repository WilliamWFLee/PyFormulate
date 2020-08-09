#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

from .models import Atom


ALIPHATIC_ORGANIC_REGEX = re.compile(r"B|C|N|O|S|P|F|Cl|Br|I")
AROMATIC_ORGANIC = "bcnosp"


class ParserError(ValueError):
    def __init__(self, msg, doc, pos):
        error_message = "{}: character {!r}, position {}".format(msg, doc[pos], pos)
        super().__init__(error_message)


class Parser:
    """
    Class for parsing SMILES
    """

    def __init__(self, formula: str):
        self.formula = formula

    def _parse_aromatic_organic(self, formula: str, start: int):
        symbol = formula[start]
        if symbol not in AROMATIC_ORGANIC:
            raise ParserError("Unknown organic aromatic element", formula, start)
        atom = Atom(symbol, aromatic=True)
        self._atom_stack.append(atom)

        return atom, start

    def _parse_aliphatic_organic(self, formula: str, start: int):
        match = ALIPHATIC_ORGANIC_REGEX.match(formula, start)
        if not match:
            raise ParserError("Unknown organic aliphatic element", formula, start)

        symbol = match.group()
        atom = Atom(symbol)
        self._atom_stack.append(atom)

        return atom, match.end() - 1

    def _parse_organic_atom(self, formula: str, start: int):
        if formula[start].islower():
            atom, end = self._parse_aromatic_organic(formula, start)
        else:
            atom, end = self._parse_aliphatic_organic(formula, start)

        return atom, end

    def parse(self):
        species_list = []
        self._atom_stack = []
        idx = 0

        species = []
        while idx < len(self.formula):
            char = self.formula[idx]
            if char in " \t\r\n":
                break
            if char.isalpha():
                atom, idx = self._parse_organic_atom(self.formula, idx)
                if self._atom_stack:
                    atom.bond(self._atom_stack[-1])
                self._atom_stack.append(atom)
                species.append(atom)

            idx += 1

        if species:
            species_list.append(species)

        del self._atom_stack
        return species_list
