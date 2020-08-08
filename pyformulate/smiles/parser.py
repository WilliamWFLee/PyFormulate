#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import Optional


class Parser:
    """
    Class for parsing SMILES
    """

    def __init__(self, formula: Optional[str] = None):
        if formula is not None:
            self.parse(formula)

    def parse(self):
        pass
