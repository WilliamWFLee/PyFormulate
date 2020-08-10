#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import List, Optional, Tuple

from .models import Atom, BondType, Molecule, Element

ALIPHATIC_ORGANIC = ("B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I")
AROMATIC_ORGANIC = "bcnosp"
AROMATIC = ("b", "c", "n", "o", "p", "s", "se", "as")
CHAR_TO_BOND_TYPE = {
    "-": BondType.SINGLE,
    "=": BondType.DOUBLE,
    "#": BondType.TRIPLE,
    "$": BondType.QUADRUPLE,
    ":": BondType.AROMATIC,
    "/": BondType.BOTTOM_TOP,
    "\\": BondType.TOP_BOTTOM,
}


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
        atom = None
        symbol = self._stream.next
        if symbol.isalpha() and symbol.islower():
            if symbol not in AROMATIC_ORGANIC:
                raise DecodeError("Unknown aromatic organic", symbol, self._stream.pos)
            atom = Atom(next(self._stream), aromatic=True)
        return atom

    def _parse_aliphatic_organic(self) -> Atom:
        symbol = self._stream.next
        if not (symbol.isalpha() and symbol.isupper()):
            return None

        next(self._stream)
        next_char = self._stream.next
        if next_char.isalpha() and next_char.islower():
            symbol += next(self._stream)
        if symbol in ALIPHATIC_ORGANIC:
            return Atom(symbol)
        raise DecodeError("Unknown aliphatic organic", symbol, self._stream.pos - 1)

    def _parse_organic_atom(self) -> Atom:
        atom = None
        for atom_parser in (
            self._parse_aromatic_organic,
            self._parse_aliphatic_organic,
        ):
            atom = atom_parser()
            if atom is not None:
                break
        return atom

    def _parse_isotope(self) -> int:
        number = ""
        while self._stream.next.isnumeric():
            number += next(self._stream)

        return int(number) if number else None

    def _parse_symbol(self) -> Tuple[Optional[Element], Optional[bool]]:
        symbol = next(self._stream)
        if symbol == "*":
            return Element.UNKNOWN, False
        if self._stream.next.isalpha() and self._stream.next.islower():
            symbol += next(self._stream)
        if symbol in AROMATIC:
            return Element[symbol.title()], True
        if symbol in Element.__members__:
            return Element[symbol], False
        return None, None

    def _parse_bracket_atom(self) -> Atom:
        open_bracket = self._stream.next
        if open_bracket != "[":
            return None

        next(self._stream)
        isotope = self._parse_isotope()
        element, aromatic = self._parse_symbol()
        if element is None:
            raise DecodeError(
                "Element symbol expected for bracket atom",
                self._stream.next,
                self._stream.pos,
            )

        atom = Atom(element, isotope=isotope, aromatic=aromatic)
        close_bracket = next(self._stream)
        if close_bracket != "]":
            raise DecodeError(
                "Unexpected character, expected ']'",
                close_bracket,
                self._stream.pos - 1,
            )

        return atom

    def _parse_atom(self) -> Optional[Atom]:
        atom = None
        for atom_parser in (self._parse_bracket_atom, self._parse_organic_atom):
            atom = atom_parser()
            if atom is not None:
                break
        return atom

    def _parse_branched_atom(self) -> Optional[Atom]:
        return self._parse_atom()

    def _parse_bond(self) -> BondType:
        bond_char = self._stream.next
        if bond_char in CHAR_TO_BOND_TYPE:
            return CHAR_TO_BOND_TYPE[next(self._stream)]
        return None

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

        if not self._stream.next:
            return atom

        if self._stream.next == ".":
            next(self._stream)
            next_atom = self._parse_chain()
            if not next_atom:
                raise DecodeError("Expected chain after dot", ".", self._stream.pos)
            return

        bond_type = self._parse_bond()
        if bond_type is not None:
            next_atom = self._parse_chain(molecule_idx)
            if next_atom is None:
                raise DecodeError(
                    "Expected atom after bond symbol",
                    self._stream.next,
                    self._stream.pos,
                )

            atom.bond(next_atom, bond_type)
            return

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
                        "Unexpected character", terminator, self._stream.pos - 1
                    )
            except StopIteration:
                pass

        return DecodeResult(
            [Molecule(molecule) for molecule in self._molecules], self._stream.remainder
        )
