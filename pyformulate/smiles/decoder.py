#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import warnings
from typing import List, Optional, Tuple

from .models import Atom, BondType, ChiralClass, Element, Molecule

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
BOND_TYPE_TO_CHAR = {v: k for k, v in CHAR_TO_BOND_TYPE.items()}
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


class DecodeWarning(Warning):
    """
    Class of all SMILES decode warnings
    """

    def __init__(self, msg, snippet, pos):
        warning_message = "{}: {!r}, position {}".format(msg, snippet, pos)
        super().__init__(warning_message)


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

    def _parse_aromatic_organic(self) -> Optional[Atom]:
        atom = None
        symbol = self._stream.next
        if symbol.isalpha() and symbol.islower():
            if symbol not in AROMATIC_ORGANIC:
                raise DecodeError("Unknown aromatic organic", symbol, self._stream.pos)
            atom = Atom(next(self._stream), aromatic=True)
        return atom

    def _parse_aliphatic_organic(self) -> Optional[Atom]:
        symbol = self._stream.next
        if not (symbol.isalpha() and symbol.isupper()):
            return None

        next(self._stream)
        next_char = self._stream.next
        if (
            next_char.isalpha()
            and next_char.islower()
            and symbol + next_char in ALIPHATIC_ORGANIC
        ):
            symbol += next(self._stream)
        if symbol in ALIPHATIC_ORGANIC:
            return Atom(symbol)
        raise DecodeError("Unknown aliphatic organic", symbol, self._stream.pos - 1)

    def _parse_organic_atom(self) -> Optional[Atom]:
        atom = None
        for atom_parser in (
            self._parse_aromatic_organic,
            self._parse_aliphatic_organic,
        ):
            atom = atom_parser()
            if atom is not None:
                break
        if atom is not None:
            atom._organic = True  # For implicit hydrogens
        return atom

    def _parse_number(self) -> Optional[int]:
        number = ""
        while self._stream.next.isnumeric():
            number += next(self._stream)
        return int(number) if number else None

    def _parse_isotope(self) -> Optional[int]:
        return self._parse_number()

    def _parse_symbol(self) -> Tuple[Optional[Element], Optional[bool]]:
        symbol = self._stream.next
        if symbol == "*":
            next(self._stream)
            return Element.UNKNOWN, False
        if not symbol.isalpha():
            return None, None
        if symbol.islower():
            if symbol in AROMATIC:
                next(self._stream)
                return Element[symbol.title()], True
            else:
                raise DecodeError("Unknown aromatic", symbol, self._stream.pos)

        next(self._stream)
        if (
            self._stream.next.islower()
            and symbol + self._stream.next in Element.__members__
        ):
            symbol += next(self._stream)
        if symbol in Element.__members__:
            return Element[symbol], False
        else:
            raise DecodeError("Unknown element", symbol, self._stream.pos - len(symbol))

    def _parse_digit(self) -> Optional[str]:
        if self._stream.next.isnumeric():
            return next(self._stream)
        return None

    def _parse_chiral(self) -> Optional[ChiralClass]:
        at_sym = self._stream.next
        if at_sym != "@":
            return None

        next(self._stream)
        if self._stream.next == "@":
            next(self._stream)
            return ChiralClass._CLOCKWISE
        if self._stream.next.isalpha() and self._stream.next.isupper():
            chiral_type = next(self._stream)
            if not (self._stream.next.isalpha() and self._stream.next.isupper()):
                raise DecodeError(
                    "Incomplete chiral class", chiral_type, self._stream.pos - 1
                )
            chiral_type += next(self._stream)

            class_number = ""
            while True:
                digit = self._parse_digit()
                if digit is None:
                    break
                class_number += digit

            if class_number == "":
                raise DecodeError(
                    "Chiral class number not specified",
                    chiral_type,
                    self._stream.pos - len(chiral_type),
                )

            chiral_type += class_number
            if chiral_type not in ChiralClass.__members__:
                raise DecodeError(
                    "Unknown chiral class",
                    chiral_type,
                    self._stream.pos - len(chiral_type),
                )
            return ChiralClass[chiral_type]

        return ChiralClass._ANTICLOCKWISE

    def _parse_hydrogen_count(self) -> int:
        if self._stream.next != "H":
            return 0
        next(self._stream)

        digit = self._parse_digit()
        if digit is None:
            return 1
        return int(digit)

    def _parse_sign(self) -> str:
        if self._stream.next not in ("+", "-"):
            return None
        return next(self._stream)

    def _parse_charge(self) -> int:
        sign = self._parse_sign()
        if sign is None:
            return 0
        next_sym = self._parse_sign()
        if next_sym is not None:
            if next_sym == sign:
                warnings.warn(
                    DecodeWarning(
                        f"Use of {2 * sign} is deprecated. Use {sign}2 instead",
                        2 * sign,
                        self._stream.pos - 2,
                    )
                )
                return 2 if sign == "+" else -2
            else:
                raise DecodeError(
                    f"Unexpected character following charge sign {sign}",
                    next_sym,
                    self._stream.pos - 1,
                )

        digit = self._parse_digit()
        if digit is None:
            return 1 if sign == "+" else -1
        next_digit = self._parse_digit()
        if next_digit is None:
            return int(digit) if sign == "+" else -int(digit)
        count = digit + next_digit
        return int(count) if sign == "+" else -int(count)

    def _parse_atom_class(self) -> int:
        if self._stream.next != ":":
            return 0
        next(self._stream)
        number = self._parse_number()
        if number is None:
            raise DecodeError(
                "Expected atom class after ':'", self._stream.next, self._stream.pos
            )
        return number

    def _parse_bracket_atom(self) -> Tuple[Optional[Atom], int]:
        open_bracket = self._stream.next
        if open_bracket != "[":
            return None, 0

        next(self._stream)
        isotope = self._parse_isotope()
        element, aromatic = self._parse_symbol()
        if element is None:
            raise DecodeError(
                "Element symbol expected for bracket atom",
                self._stream.next,
                self._stream.pos,
            )
        chiral_class = self._parse_chiral()
        hydrogen_count = self._parse_hydrogen_count()
        if hydrogen_count > 0 and element == Element.H:
            raise DecodeError(
                "Hydrogen count for bracket atom of hydrogen not permitted",
                "H",
                self._stream.pos - 2,
            )
        charge = self._parse_charge()
        atom_class = self._parse_atom_class()

        atom = Atom(
            element,
            isotope=isotope,
            aromatic=aromatic,
            chiral_class=chiral_class,
            charge=charge,
            atom_class=atom_class,
        )

        close_bracket = next(self._stream)
        if close_bracket != "]":
            raise DecodeError(
                "Unexpected character, expected ']'",
                close_bracket,
                self._stream.pos - 1,
            )

        return atom, hydrogen_count

    def _parse_atom(self) -> Tuple[Optional[Atom], int]:
        atom, hydrogen_count = self._parse_bracket_atom()
        if atom:
            return atom, hydrogen_count

        organic_atom = self._parse_organic_atom()
        if organic_atom is not None:
            return organic_atom, 0
        return None, 0

    def _parse_bond(self) -> Optional[BondType]:
        bond_char = self._stream.next
        if bond_char in CHAR_TO_BOND_TYPE:
            return CHAR_TO_BOND_TYPE[next(self._stream)]
        return None

    def _parse_ring_bond(self) -> Tuple[BondType, int]:
        bond_type = self._parse_bond()
        digit = self._parse_digit()
        if digit is not None:
            return bond_type, int(digit)
        if self._stream.next == "%":
            next(self._stream)
            digit = self._parse_digit()
            if digit is None:
                raise DecodeError(
                    "Expected digit after '%' for ring number",
                    self._stream.next,
                    self._stream.pos,
                )
            next_digit = self._parse_digit()
            if next_digit is None:
                raise DecodeError(
                    "Expected two-digit number after '%'",
                    self._stream.next,
                    self._stream.pos,
                )
            return bond_type, int(digit + next_digit)
        return bond_type, None

    def _parse_branch(
        self, molecule: Molecule
    ) -> Optional[Tuple[Optional[Atom], Optional[BondType]]]:
        if self._stream.next != "(":
            return None
        next(self._stream)

        atom = None
        bond_type = self._parse_bond()
        if self._stream.next != ".":
            if bond_type is None:
                bond_type = BondType.SINGLE
            atom = self._parse_chain(molecule)
        else:
            next(self._stream)
            atom = self._parse_chain()
        if atom is None:
            raise DecodeError(
                "Expected chain within branch", self._stream.next, self._stream.pos
            )
        if self._stream.next != ")":
            raise DecodeError(
                "Expected ')' to close branch", self._stream.next, self._stream.pos
            )
        next(self._stream)
        return atom, bond_type

    def _parse_branched_atom(self, molecule: Molecule) -> Tuple[Atom, int, BondType]:
        atom, hydrogen_count = self._parse_atom()
        while True:
            bond_type, rnum = self._parse_ring_bond()
            if rnum is None:
                break
            if rnum in self._rnums:
                other_atom, other_bond_type = self._rnums[rnum]
                if (
                    bond_type is not None
                    and other_bond_type is not None
                    and bond_type != other_bond_type
                ):
                    raise DecodeError(
                        (
                            f"Mismatched bond type for rnum {rnum}, "
                            f"expected {BOND_TYPE_TO_CHAR[other_bond_type]}"
                        ),
                        BOND_TYPE_TO_CHAR[bond_type],
                        self._stream.pos,
                    )
                bond_type = bond_type if bond_type is not None else other_bond_type
                atom.bond(other_atom, bond_type)
                del self._rnums[rnum]
            else:
                self._rnums[rnum] = (atom, bond_type)
        branches = False
        bond_pos = self._stream.pos - 1
        while True:
            result = self._parse_branch(molecule)
            if result is None:
                break
            if bond_type is not None:
                raise DecodeError(
                    "Unexpected bond symbol before branch",
                    BOND_TYPE_TO_CHAR[bond_type],
                    bond_pos,
                )
            branches = True
            next_atom, next_atom_bond_type = result
            if next_atom is not None and next_atom_bond_type is not None:
                atom.bond(next_atom, next_atom_bond_type)

        if branches:
            bond_type = self._parse_bond()
        return atom, hydrogen_count, bond_type

    def _parse_chain(self, molecule: Optional[Molecule] = None) -> Atom:
        if molecule is None:
            molecule = Molecule()
            self._molecules.append(molecule)

        atom, hydrogen_count, bond_type = self._parse_branched_atom(molecule)
        if atom:
            molecule.add(atom)
            for _ in range(hydrogen_count):
                molecule.new_bonded_atom(atom, bond_type, Element.H)
        else:
            return None

        if self._stream.next == ".":
            next(self._stream)
            next_atom = self._parse_chain()
            if not next_atom:
                raise DecodeError("Expected chain after dot", ".", self._stream.pos)
            return None
        if bond_type is not None:
            next_atom = self._parse_chain(molecule)
            if next_atom is None:
                raise DecodeError(
                    "Expected atom after bond symbol",
                    self._stream.next,
                    self._stream.pos,
                )
            molecule.bond(atom, next_atom, bond_type)
        else:
            next_atom = self._parse_chain(molecule)
            if next_atom is not None:
                molecule.bond(atom, next_atom)
        return atom

    def _add_implicit_hydrogens(self, molecule: Molecule):
        for atom in molecule:
            if (
                atom.element in VALENCIES
                and getattr(atom, "_organic", False)
                and atom.valency <= max(VALENCIES[atom.element])
            ):
                for valency in VALENCIES[atom.element]:
                    valency_diff = valency - atom.valency
                    if valency_diff > 0:
                        for _ in range(valency_diff):
                            molecule.new_bonded_atom(atom, element=Element.H)
                        break
                    elif valency_diff == 0:
                        break
                del atom._organic

    def decode(self) -> DecodeResult:
        self._molecules = []
        self._rnums = {}
        self._stream = PeekableStream(self.smiles)
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
            if self._rnums:
                raise DecodeError(
                    (
                        "Unexpected end-of-input with unmatched rnums "
                        + ", ".join(str(k) for k in self._rnums.keys())
                    ),
                    self._stream.next,
                    self._stream.pos,
                )

        for molecule in self._molecules:
            self._add_implicit_hydrogens(molecule)
        try:
            return DecodeResult(self._molecules, self._stream.remainder)
        finally:
            del self._molecules


def loads(s: str) -> DecodeResult:
    """
    Parses a SMILES string, and produces a :class:`DecodeResult` object.

    :param s: The SMILES string
    :type s: str
    :return: The result of decoding the SMILES
    :rtype: DecodeResult
    """
    return Decoder(s).decode()
