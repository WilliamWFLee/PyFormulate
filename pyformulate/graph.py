#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pyformulate.graph

Module for graph-based functionality
"""

from collections import defaultdict
from typing import Hashable, Any


class InfoGraph:
    """
    A class for defining a graph whose edges have additional information
    associated with them.

    Programmatically, the graph is implemented as a dictionary of dictionaries.
    The keys of the outer dictionary are the nodes of the graphs,
    and the keys of each inner dictionary are the nodes to which to the corresponding
    node is attached to. The values of the inner dictionaries
    is the value associated with that particular edge.

    For example, the dictionary:

    .. code-block:: py

        {
            "A": {"B": 1},
            "B": {"A": 1},
        }

    represents the graph with nodes `"A"` and `"B"`, with a single edge between them,
    having the associated value of `1`, representing a weight.
    """

    def __init__(self):
        self._dict = defaultdict(dict)

    def connect(self, value: Hashable, other: Hashable, info: Any = None):
        """
        Connects two values representing nodes on the graph together,
        with the optional info describing the edge

        :param value: One of the node's value
        :type value: Hashable
        :param other: The other node's value
        :type other: Hashable
        :param info: The info associated with the edge between them, defaults to None
        :type info: Any
        """
        self._dict[value][other] = info
        self._dict[other][value] = info

    def is_connected(self, value, other):
        """
        Determines whether or not two nodes are connected on the graph

        :param value: The value of one node
        :type value: Hashable
        :param other: The value of the other node
        :type other: Hashable
        :return: Whether the nodes are connected
        :rtype: bool
        """
        try:
            self._dict[value][other]
            return True
        except KeyError:
            return False

    def __contains__(self, value):
        return value in self._dict
