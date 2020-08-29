#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import warnings
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

from ..models import BondType, Element
from .models import (
    ALIPHATIC_ORGANIC,
    AROMATIC,
    AROMATIC_ORGANIC,
    CHAR_TO_BOND_ORDER,
    CHAR_TO_CIS_TRANS,
    VALENCIES,
    Atom,
    ChiralClass,
    Molecule,
)


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
        """
        The next character in the stream
        """
        try:
            return self.value[self.pos]
        except IndexError:
            return ""

    @property
    def remainder(self):
        """
        The remainder of the string in the stream
        """
        return self.value[self.pos :]

    def __next__(self):
        try:
            return self.value[self.pos]
        except IndexError:
            # End of string reached, raises StopIteration as defined
            raise StopIteration() from None
        finally:
            self.pos += 1

    def __iter__(self):
        return self


class DecodeError(ValueError):
    """
    Class of all SMILES decode errors

    Inherits from :class:`ValueError`
    """

    def __init__(
        self, msg: str, snippet: Optional[str] = None, pos: Optional[int] = None
    ):
        """
        Creates an instance of a :class:`DecodeError`.

        Has three attributes: the error's message, the character/phrase of the error,
        and the position the error occurred

        :param msg: The error message
        :type msg: str
        :param snippet: The character/phrase causing the error, defaults to None
        :type snippet: Optional[str]
        :param pos: The position of the characters causing the error in the string, defaults to None
        :type pos: Optional[int
        """
        error_message = msg + (
            f": {snippet!r}" + (f", position {pos}" if pos is not None else "")
            if snippet is not None
            else ""
        )
        super().__init__(error_message)


class DecodeWarning(Warning):
    """
    Class of all SMILES decode warnings

    Inherits from :class:`Warning`
    """

    def __init__(self, msg, snippet, pos):
        """
        Creates an instance of a :class:`DecodeWarning`.

        Has three attributes: the warnings's message,
        the character/phrase of the warning, and the position the warning occurred

        :param msg: The warning message
        :type msg: str
        :param snippet: The character/phrase causing the warning
        :type snippet: str
        :param pos: The position of the characters causing the warning in the string
        :type pos: int
        """
        warning_message = "{}: {!r}, position {}".format(msg, snippet, pos)
        super().__init__(warning_message)


class DecodeResult:
    """
    Represents the result of decoding a SMILES string

    .. attribute:: molecules

        The molecules produced from the SMILES string

        :type: List[Molecule]

    .. attribute:: remainder

        The rest of the SMILES string after the terminator of the SMILES

        :type: str
    """

    def __init__(self, molecules: List[Molecule], remainder: str):
        self.molecules = molecules
        self.remainder = remainder


class Decoder:
    """
    Class for decoding SMILES

    .. attribute:: smiles

        The SMILES string being decoded

        :type: str
    """

    def __init__(self, smiles: str):
        """
        Create a new instance of :class:`Decoder` to parse a given SMILES string

        :param smiles: The SMILES string to parse
        :type smiles: str
        """
        self.smiles = smiles

    def _parse_aromatic_organic(self) -> Optional[Atom]:
        # Aromatic organics are single-character only
        atom = None
        symbol = self._stream.next
        if symbol.isalpha() and symbol.islower():
            if symbol not in AROMATIC_ORGANIC:
                raise DecodeError("Unknown aromatic organic", symbol, self._stream.pos)
            atom = Atom(element=next(self._stream), aromatic=True)
        return atom

    def _parse_aliphatic_organic(self) -> Optional[Atom]:
        symbol = self._stream.next
        if not (symbol.isalpha() and symbol.isupper()):
            return None

        next(self._stream)
        next_char = self._stream.next
        if (  # If the next character forms a valid two-character aliphatic organic
            next_char.isalpha()
            and next_char.islower()
            and symbol + next_char in ALIPHATIC_ORGANIC
        ):
            symbol += next(self._stream)
        if symbol in ALIPHATIC_ORGANIC:
            return Atom(element=symbol)
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
        if self._stream.next == "@":  # @@
            next(self._stream)
            return ChiralClass._CLOCKWISE
        if self._stream.next.isalpha() and self._stream.next.isupper():
            chiral_type = next(self._stream)
            # Chiral specification are two uppercase letters
            if not (self._stream.next.isalpha() and self._stream.next.isupper()):
                raise DecodeError(
                    "Incomplete chiral class", chiral_type, self._stream.pos - 1
                )
            chiral_type += next(self._stream)

            # Chiral class number is interpreted as an arbitrary-length integer
            # TODO: Change to accept one or two digits
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
                # Use of ++ and -- is deprecated
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
            return 1 if sign == "+" else -1  # Implied magnitude of 1
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
            return None, 0  # Not a bracket atom

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
            # Hydrogen count for hydrogen is not permitted
            # under the SMILES specification
            raise DecodeError(
                "Hydrogen count for bracket atom of hydrogen not permitted",
                "H",
                self._stream.pos - 2,
            )
        charge = self._parse_charge()
        atom_class = self._parse_atom_class()

        atom = Atom(
            element=element,
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
            return organic_atom, 0  # Hydrogen count is implied
        return None, 0  # No atom

    def _parse_bond(self) -> Optional[BondType]:
        bond_char = self._stream.next
        if bond_char not in CHAR_TO_BOND_ORDER:
            return None
        bond_order = CHAR_TO_BOND_ORDER[next(self._stream)]
        aromatic = bond_char == ":"
        cis_trans = CHAR_TO_CIS_TRANS.get(bond_char, None)

        return BondType(bond_order, aromatic, cis_trans)

    def _parse_ring_bond(self) -> Tuple[Optional[BondType], int]:
        bond_type = self._parse_bond()
        digit = self._parse_digit()
        if digit is not None:
            return bond_type, int(digit)
        if self._stream.next == "%":  # Two-digit number
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
        return bond_type, None  # Standard bond

    def _parse_branch(
        self, molecule: Molecule
    ) -> Optional[Tuple[Optional[Atom], Optional[BondType]]]:
        if self._stream.next != "(":
            return None  # Not the start of a branch
        next(self._stream)

        atom = None
        bond_type = self._parse_bond()
        if self._stream.next != ".":  # Standard branch
            atom = self._parse_chain(molecule)
        else:
            next(self._stream)
            atom = self._parse_chain()  # Disconnected molecule branch
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

    def _parse_branched_atom(
        self, molecule: Molecule
    ) -> Tuple[Optional[Atom], int, Optional[BondType]]:
        atom, hydrogen_count = self._parse_atom()
        if atom is None:
            return None, 0, None
        molecule.add(atom)
        # Adds the hydrogens
        for _ in range(hydrogen_count):
            hydrogen = molecule.new_atom(element=Element.H)
            molecule.bond(atom, hydrogen)
        # Parses ring bond and standard bond
        while True:
            bond_type, rnum = self._parse_ring_bond()
            if rnum is None:  # Not a ring bond
                break
            if rnum in self._rnums:
                # Completes the ring closure
                other_atom, other_bond_type = self._rnums[rnum]
                if (
                    bond_type is not None
                    and other_bond_type is not None
                    and bond_type != other_bond_type
                ):  # Ring bond types must match if both have specified it
                    raise DecodeError(
                        f"Mismatched bond type for rnum {rnum}", "", self._stream.pos,
                    )
                bond_type = bond_type if bond_type is not None else other_bond_type
                if other_atom.molecule != molecule:
                    self._molecules.remove(other_atom.molecule)
                    molecule.merge(other_atom.molecule)
                atom.bond(other_atom, bond_type)
                del self._rnums[rnum]  # Delete rnum for reuse
            else:  # Rnum not already encountered
                self._rnums[rnum] = (atom, bond_type)
        branches = False  # If this atom has branches
        bond_pos = self._stream.pos - 1  # For error-raising
        while True:
            result = self._parse_branch(molecule)
            if result is None:
                break  # No more branches
            if bond_type is not None:  # If bond symbol is encountered before branch
                raise DecodeError(
                    "Unexpected bond symbol before branch",
                    self._stream.value[bond_pos],
                    bond_pos,
                )
            branches = True
            next_atom, next_atom_bond_type = result
            if next_atom is not None and next_atom_bond_type is not None:
                atom.bond(next_atom, next_atom_bond_type)  # Bonds the atoms

        if branches:
            bond_type = self._parse_bond()  # Parse the bond symbol after the branches
        return atom, hydrogen_count, bond_type

    def _parse_chain(self, molecule: Optional[Molecule] = None) -> Atom:
        # Creates a new molecule if not specified
        if molecule is None:
            molecule = Molecule()
            self._molecules.append(molecule)

        # Parse an atom and its branches
        atom, hydrogen_count, bond_type = self._parse_branched_atom(molecule)
        if atom is None:
            return None  # Signals end of chain

        if self._stream.next == ".":  # Disconnected molecule
            next(self._stream)
            next_atom = self._parse_chain()
            if not next_atom:
                raise DecodeError("Expected chain after dot", ".", self._stream.pos)
        elif bond_type is not None:  # An atom is expected afterwards
            next_atom = self._parse_chain(molecule)
            if next_atom is None:
                raise DecodeError(
                    "Expected atom after bond symbol",
                    self._stream.next,
                    self._stream.pos,
                )
            molecule.bond(atom, next_atom, bond_type)
        else:  # Implicit single bond, or no next atom
            next_atom = self._parse_chain(molecule)
            if next_atom is not None:
                molecule.bond(atom, next_atom)
        return atom

    def _add_implicit_hydrogens(self, molecule: Molecule):
        for atom in molecule:
            # Atom must be organic, have its element have defined "normal" valencies,
            # and its current valency must be less than the maximum "normal" valency
            # for its element
            if (
                atom.element in VALENCIES
                and getattr(atom, "_organic", False)
                and atom.valency <= max(VALENCIES[atom.element])
            ):
                for valency in VALENCIES[atom.element]:
                    # Aromatic atoms receive one fewer hydrogen
                    valency_diff = valency - atom.valency - (1 if atom.aromatic else 0)
                    if valency_diff > 0:
                        # If the current valency is lower than one of the valencies
                        for _ in range(valency_diff):
                            hydrogen = molecule.new_atom(element=Element.H)
                            molecule.bond(atom, hydrogen)
                        break  # Ensures valency is always next-highest
                    elif valency_diff == 0:
                        break  # If the atom is already one of the valencies
                del atom._organic  # Cleans the attribute from the atom

    def _determine_double_bonds(
        self,
        bonds: Dict[Atom, Set[Atom]],
        double_bonds: Optional[Dict[Atom, Set[Atom]]] = None,
    ) -> Optional[Dict[Atom, Set[Atom]]]:
        # Finds a perfect matching for the graph given by the atoms and bonds
        if not bonds:
            return double_bonds
        if double_bonds is None:
            double_bonds = defaultdict(set)
        for atom, neighbours in bonds.items():
            for other_atom in neighbours:
                proposed_double_bonds = self._determine_double_bonds(
                    {
                        k: {p for p in v if p not in (atom, other_atom)}
                        for k, v in bonds.items()
                        if k not in (atom, other_atom)
                    },
                    double_bonds,
                )
                if proposed_double_bonds is not None:
                    proposed_double_bonds[atom].add(other_atom)
                    return proposed_double_bonds
        return None

    def _kekulize(self, molecule: Molecule):
        def aromatic_predicate(x):
            return x.aromatic and (
                x.valency not in VALENCIES[x.element]
                if x.element in VALENCIES
                else True
            )

        bonds = {
            k: {p for p in v if aromatic_predicate(p)}
            for k, v in molecule.bonds.items()
            if aromatic_predicate(k)
        }
        if not bonds:
            return
        double_bonds = self._determine_double_bonds(bonds)
        if double_bonds is None:
            raise DecodeError(
                f"Failed to verify aromaticity of molecule with formula {molecule}"
            )

    def decode(self) -> DecodeResult:
        """
        Decode the SMILES string of this decoder.

        :raises DecodeError: If the SMILES is syntactically or semantically incorrect
        :return: The result of decoding the SMILES string
        :rtype: DecodeResult
        """
        self._molecules = []
        self._rnums = {}
        self._stream = PeekableStream(self.smiles)

        # If the SMILES doesn't start with a terminator
        if self._stream.next not in " \t\r\n":
            try:
                # Parses the SMILES string
                self._parse_chain()
            except StopIteration:
                raise DecodeError("Unexpected end-of-input", "", self._stream.pos)
            try:
                terminator = next(self._stream)
                # Expects a terminator
                if terminator not in " \t\r\n":
                    raise DecodeError(
                        "Unexpected character", terminator, self._stream.pos - 1
                    )
            except StopIteration:
                # End-of-string is accepted as a terminator
                pass
            if self._rnums:
                # Ring bond numbers must be matched
                raise DecodeError(
                    (
                        "Unexpected end-of-input with unmatched rnums "
                        + ", ".join(str(k) for k in self._rnums.keys())
                    ),
                    self._stream.next,
                    self._stream.pos,
                )

        # Adds implicit hydrogens to organic atoms
        for molecule in self._molecules:
            self._add_implicit_hydrogens(molecule)
            self._kekulize(molecule)
        try:
            return DecodeResult(self._molecules, self._stream.remainder)
        finally:
            # Clean up namespace
            del self._molecules
            del self._stream


def loads(s: str) -> DecodeResult:
    """
    Parses a SMILES string, and produces a :class:`DecodeResult` object.

    :param s: The SMILES string
    :type s: str
    :return: The result of decoding the SMILES
    :rtype: DecodeResult
    """
    return Decoder(s).decode()
