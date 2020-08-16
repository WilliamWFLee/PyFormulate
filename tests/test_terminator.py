#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pyformulate.smiles import loads


def test_terminator():
    for terminator in " \t\r\n":
        assert not loads(terminator).molecules
