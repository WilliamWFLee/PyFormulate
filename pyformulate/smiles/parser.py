#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import Optional
from enum import Enum, auto


class Symbol(Enum):
    TERMINATOR = auto()
    UCASE_ALPHA = auto()
    LCASE_ALPHA = auto()
    DIGIT = auto()
    OPEN_PAREN = auto()
    CLOSE_PAREN = auto()
    OPEN_BRACK = auto()
    CLOSE_BRACK = auto()
    AT = auto()
    DOT = auto()
    HYPHEN = auto()
    PLUS = auto()
    EQUALS = auto()
    HASH = auto()
    DOLLAR = auto()
    COLON = auto()
    F_SLASH = auto()
    B_SLASH = auto()
    PERCENT = auto()
    ASTERISK = auto()
    UNKNOWN = auto()


class Parser:
    """
    Class for parsing SMILES
    """

    def __init__(self, formula: Optional[str] = None):
        if formula is not None:
            self.parse(formula)

    def parse(self):
        pass
