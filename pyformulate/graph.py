#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pyformulate.graph

Module for graph-based functionality
"""

from collections import defaultdict


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
        self.graph = defaultdict(dict)

    def connect(self, value, other, info=None):
        """
        Connects two values representing nodes on the graph together,
        with the optional info describing the edge

        :param value: One of the node's value
        :type value: Any
        :param other: The other node's value
        :type other: Any
        :param info: The info associated with the edge between them, defaults to None
        :type info: Any
        """
        self.graph[value][other] = info
        self.graph[other][value] = info
