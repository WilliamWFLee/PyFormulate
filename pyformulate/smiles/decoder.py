#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from typing import List, Optional, Tuple

from .models import Atom, BondType, Molecule

ALIPHATIC_ORGANIC_REGEX = re.compile(r"B|C|N|O|S|P|F|Cl|Br|I")
AROMATIC_ORGANIC = "bcnosp"


class PeekableStream:
    """
    Class for producing a peekable stream of characters from a string,
    that can peek one character ahead.
    """

    def __init__(self, value: str):
        self.value = value
        self.pos = 0
        self._fill()

    @property
    def next(self):
        try:
            return self.value[self.pos]
        except IndexError:
            return None

    def __next__(self):
        try:
            return self.value[self.pos]
            self.pos += 1
        except IndexError:
            raise StopIteration() from None

    def __iter__(self):
        return self


class DecodeError(ValueError):
    def __init__(self, msg, doc, pos):
        error_message = "{}: character {!r}, position {}".format(msg, doc[pos], pos)
        super().__init__(error_message)


class DecodeResult:
    """
    Represents the result of decoding a SMILES string
    """

    def __init__(self, molecules: List[Molecule], remainder: str):
        self.molecules = molecules
        self.remainder = remainder


class Decoder:
    """
    Class for decoding SMILES
    """

    def __init__(self, smiles: str):
        self.smiles = smiles
        self.species_list = []

    def _parse_aromatic_organic(self, start: int) -> Tuple[Atom, int]:
        symbol = self.smiles[start]
        if symbol not in AROMATIC_ORGANIC:
            raise DecodeError("Unknown organic aromatic element", self.smiles, start)
        atom = Atom(symbol, aromatic=True)
        return atom, start

    def _parse_aliphatic_organic(self, start: int) -> Tuple[Atom, int]:
        match = ALIPHATIC_ORGANIC_REGEX.match(self.smiles, start)
        if not match:
            raise DecodeError("Unknown organic aliphatic element", self.smiles, start)

        symbol = match.group()
        atom = Atom(symbol)
        return atom, match.end() - 1

    def _parse_organic_atom(self, start: int) -> Tuple[Atom, int]:
        if self.smiles[start].islower():
            atom, end = self._parse_aromatic_organic(start)
        else:
            atom, end = self._parse_aliphatic_organic(start)

        return atom, end

    def _parse_branched_atom(self, start: int) -> Tuple[Optional[Atom], int]:
        char = self.smiles[start : start + 1]
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
            return char_to_bond_type[self.smiles[start]], start
        except KeyError:
            raise DecodeError("Unknown bond type", self.smiles, start)

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
            molecule.append(atom)
            chain.append(atom)
            idx += 1
        else:
            return molecule_idx, idx

        char = self.smiles[idx : idx + 1]
        if char:
            if char in r"\/-=#$:":
                bond_type, idx = self._parse_bond(idx)
                idx += 1

                length_before = len(chain)
                molecule_idx, idx = self._parse_chain(idx, molecule_idx, chain)
                length_after = len(chain)

                if length_before == length_after:
                    raise DecodeError(
                        "Expected atom after bond symbol", self.smiles, idx - 1
                    )

                last_atom = chain.pop()
                chain[-1].bond(last_atom, bond_type)
            elif char == ".":
                idx += 1
                molecule_idx, idx = self._parse_chain(idx)
            else:
                length_before = len(chain)
                molecule_idx, idx = self._parse_chain(idx, molecule_idx, chain)
                length_after = len(chain)

                if length_before != length_after:
                    last_atom = chain.pop()
                    chain[-1].bond(last_atom)
                    idx += 1

        return molecule_idx, idx

    def decode(self) -> Tuple[List[List[Atom]], str]:
        self._molecules = []
        idx = 0

        if self.smiles[0] not in " \t\r\n":
            _, idx = self._parse_chain(0)
            char = self.smiles[idx : idx + 1]
            if char not in " \t\r\n":
                raise DecodeError("Unexpected character", self.smiles, idx)

        remainder = self.smiles[idx + 2 :]
        return DecodeResult(
            [Molecule(molecule) for molecule in self._molecules], remainder
        )
