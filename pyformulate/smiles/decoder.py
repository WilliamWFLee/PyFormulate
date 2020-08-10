#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import List, Optional, Tuple

from .models import Atom, BondType, Molecule

ALIPHATIC_ORGANIC = ("B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I")
AROMATIC_ORGANIC = "bcnosp"


class PeekableStream:
    """
    Class for producing a peekable stream of characters from a string,
    that can peek one character ahead.
    """

    def __init__(self, value: str):
        self.value = value
        self.pos = 0

    @property
    def next(self):
        try:
            return self.value[self.pos]
        except IndexError:
            return ""

    @property
    def remainder(self):
        return self.value[self.pos :]

    def __next__(self):
        try:
            return self.value[self.pos]
        except IndexError:
            raise StopIteration() from None
        finally:
            self.pos += 1

    def __iter__(self):
        return self


class DecodeError(ValueError):
    def __init__(self, msg, snippet, pos):
        error_message = "{}: {!r}, position {}".format(msg, snippet, pos)
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
        self._stream = PeekableStream(smiles)

    def _parse_aromatic_organic(self) -> Atom:
        symbol = next(self._stream)
        if symbol not in AROMATIC_ORGANIC:
            raise DecodeError(
                "Unknown organic aromatic element", symbol, self._stream.pos - 1
            )
        atom = Atom(symbol, aromatic=True)
        return atom

    def _parse_aliphatic_organic(self) -> Atom:
        symbol = next(self._stream)
        if symbol in ALIPHATIC_ORGANIC:
            return Atom(symbol)

        next_char = self._stream.next
        if not (next_char.isalpha() and next_char.islower()):
            raise DecodeError(
                "Unknown organic aliphatic element", symbol, self._stream.pos - 1
            )

        symbol += next(self._stream)
        if symbol not in ALIPHATIC_ORGANIC:
            raise DecodeError(
                "Unknown organic aliphatic element", symbol, self._stream.pos - 2
            )
        return Atom(symbol)

    def _parse_organic_atom(self) -> Atom:
        if self._stream.next.islower():
            atom = self._parse_aromatic_organic()
        else:
            atom = self._parse_aliphatic_organic()

        return atom

    def _parse_branched_atom(self) -> Optional[Atom]:
        if self._stream.next.isalpha():
            atom = self._parse_organic_atom()
            return atom
        return None

    def _parse_bond(self) -> Tuple[BondType]:
        char_to_bond_type = {
            "-": BondType.SINGLE,
            "=": BondType.DOUBLE,
            "#": BondType.TRIPLE,
            "$": BondType.QUADRUPLE,
            ":": BondType.AROMATIC,
            "/": BondType.BOTTOM_TOP,
            "\\": BondType.TOP_BOTTOM,
        }

        bond = next(self._stream)
        try:
            return char_to_bond_type[bond]
        except KeyError:
            raise DecodeError("Unknown bond type", bond, self._stream.pos - 1)

    def _parse_chain(self, molecule_idx: Optional[int] = None) -> Atom:
        if molecule_idx is None:
            self._molecules.append([])
            molecule_idx = len(self._molecules) - 1

        molecule = self._molecules[molecule_idx]
        atom = self._parse_branched_atom()
        if atom:
            molecule.append(atom)
        else:
            return None

        char = self._stream.next
        if not char:
            return atom
        if char in r"\/-=#$:":
            bond_type = self._parse_bond()
            next_atom = self._parse_chain(molecule_idx)

            if next_atom is None:
                raise DecodeError(
                    "Expected atom after bond symbol", char, self._stream.pos
                )

            atom.bond(next_atom, bond_type)
        elif char == ".":
            next(self._stream)
            self._parse_chain()
        else:
            next_atom = self._parse_chain(molecule_idx)
            if next_atom is not None:
                atom.bond(next_atom)

    def decode(self) -> DecodeResult:
        self._molecules = []
        if self._stream.next not in " \t\r\n":
            try:
                self._parse_chain()
            except StopIteration:
                raise DecodeError("Unexpected end-of-input", "", self._stream.pos)
            try:
                terminator = next(self._stream)
                if terminator not in " \t\r\n":
                    raise DecodeError(
                        "Unexpected character", terminator, self._stream.pos
                    )
            except StopIteration:
                pass

        return DecodeResult(
            [Molecule(molecule) for molecule in self._molecules], self._stream.remainder
        )
