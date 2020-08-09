#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from enum import Enum, auto


class TokenName(Enum):
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


class Lexer:
    """
    Lexer class for SMILES
    """

    pass
