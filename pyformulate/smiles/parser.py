#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from typing import List, Optional, Tuple

from .models import Atom, BondType

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

    def _parse_aromatic_organic(self, start: int) -> Tuple[Atom, int]:
        symbol = self.formula[start]
        if symbol not in AROMATIC_ORGANIC:
            raise ParserError("Unknown organic aromatic element", self.formula, start)
        atom = Atom(symbol, aromatic=True)
        return atom, start

    def _parse_aliphatic_organic(self, start: int) -> Tuple[Atom, int]:
        match = ALIPHATIC_ORGANIC_REGEX.match(self.formula, start)
        if not match:
            raise ParserError("Unknown organic aliphatic element", self.formula, start)

        symbol = match.group()
        atom = Atom(symbol)
        return atom, match.end() - 1

    def _parse_organic_atom(self, start: int) -> Tuple[Atom, int]:
        if self.formula[start].islower():
            atom, end = self._parse_aromatic_organic(start)
        else:
            atom, end = self._parse_aliphatic_organic(start)

        return atom, end

    def _parse_branched_atom(self, start: int) -> Tuple[Optional[Atom], int]:
        char = self.formula[start : start + 1]
        if char.isalpha():
            atom, end = self._parse_organic_atom(start)
            return atom, end
        return None, start

    def _parse_bond(self, start: int) -> Tuple[BondType, int]:
        char_to_bond_type = {
            "-": BondType.SINGLE,
            "=": BondType.DOUBLE,
            "#": BondType.TRIPLE,
            "$": BondType.QUADRUPLE,
            ":": BondType.AROMATIC,
            "/": BondType.BOTTOM_TOP,
            "\\": BondType.TOP_BOTTOM,
        }

        try:
            return char_to_bond_type[self.formula[start]], start
        except KeyError:
            raise ParserError("Unknown bond type", self.formula, start)

    def _parse_chain(
        self,
        start: int = 0,
        molecule_idx: Optional[int] = None,
        chain: Optional[List[Atom]] = None,
    ) -> Tuple[int, int]:
        if chain is None:
            chain = []
        if molecule_idx is None:
            self._molecules.append([])
            molecule_idx = len(self._molecules) - 1

        molecule = self._molecules[molecule_idx]

        atom, idx = self._parse_branched_atom(start)
        if atom:
            self._molecules[molecule_idx].append(atom)
            chain.append(atom)
            idx += 1
        else:
            return molecule_idx, idx

        char = self.formula[idx : idx + 1]
        if char:
            if char in r"\/-=#$:":
                bond_type, idx = self._parse_bond(idx)
                idx += 1

                length_before = len(molecule)
                molecule_idx, idx = self._parse_chain(idx, molecule_idx, chain)
                length_after = len(molecule)

                if length_before == length_after:
                    raise ParserError(
                        "Expected atom after bond symbol", self.formula, idx - 1
                    )

                chain.pop()
                chain[-1].bond(molecule[length_before], bond_type, chain)
            else:
                atom, idx = self._parse_chain(idx, molecule_idx)
                if atom:
                    idx += 1
                    chain.append(atom)
                    chain[-1].bond(molecule[-1])

        return molecule_idx, idx

    def parse(self) -> Tuple[List[List[Atom]], str]:
        self._molecules = []
        idx = 0

        if self.formula[0] not in " \t\r\n":
            _, idx = self._parse_chain(0)
            char = self.formula[idx : idx + 1]
            if char not in " \t\r\n":
                raise ParserError("Unexpected character", self.formula, idx)

        remainder = self.formula[idx + 2 :]
        return self._molecules, remainder
