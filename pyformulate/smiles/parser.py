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
        self.species_list = []

    def _parse_aromatic_organic(self, start: int):
        symbol = self.formula[start]
        if symbol not in AROMATIC_ORGANIC:
            raise ParserError("Unknown organic aromatic element", self.formula, start)
        atom = Atom(symbol, aromatic=True)
        return atom, start

    def _parse_aliphatic_organic(self, start: int):
        match = ALIPHATIC_ORGANIC_REGEX.match(self.formula, start)
        if not match:
            raise ParserError("Unknown organic aliphatic element", self.formula, start)

        symbol = match.group()
        atom = Atom(symbol)
        return atom, match.end() - 1

    def _parse_organic_atom(self, start: int):
        if self.formula[start].islower():
            atom, end = self._parse_aromatic_organic(start)
        else:
            atom, end = self._parse_aliphatic_organic(start)

        return atom, end

    def _parse_chain(self, start: int = 0):
        atom_stack = []
        chain = []
        idx = start

        while idx < len(self.formula):
            char = self.formula[idx]
            if char.isalpha():
                atom, idx = self._parse_organic_atom(idx)
                if atom_stack:
                    atom.bond(atom_stack[-1])
                atom_stack.append(atom)
                chain.append(atom)
            else:
                return False, idx
            idx += 1

        self.species_list.append(chain)
        return True, idx

    def parse(self):
        idx = 0

        while idx < len(self.formula):
            char = self.formula[idx]
            if char in " \t\r\n":
                break
            result, idx = self._parse_chain(idx)
            if not result:
                raise ParserError("Unexpected character", self.formula, idx)

            idx += 1

        return self.species_list
