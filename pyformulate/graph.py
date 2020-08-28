#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pyformulate.graph

Module for graph-based functionality
"""

from collections import defaultdict
from typing import Any, Dict, Hashable, Optional


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

    def add(self, value: Hashable):
        """
        Adds a value as a node to this graph

        :param value: The value to add
        :type value: Hashable
        """
        self._dict[value] = {}

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

    def are_connected(self, value: Hashable, other: Hashable) -> bool:
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

    def neighbours_of(self, value: Hashable) -> Dict[Hashable, Any]:
        try:
            return self._dict[value]
        except KeyError:
            raise ValueError(f"{value!r} is not in this graph")

    def __contains__(self, value):
        return value in self._dict


class Node:
    """
    A class for representing a node on a graph.

    Use of this class to represent the nodes
    in an instance of :class:`InfoGraph` is completely optional,
    however it can be useful to derive from this class,

    .. attribute:: graph

        The instance of :class:`InfoGraph` this node is associated with
    """

    def __init__(self, graph: Optional[InfoGraph] = None):
        """
        Instantiates an instance of a node, with an optional parameter
        for the instance :class:`InfoGraph` it is associated to.
        If it is not specified, a new graph is created for it.

        :param graph: The graph, defaults to None
        :type graph: Optional[InfoGraph]
        """
        if graph is None:
            graph = InfoGraph()
        self.graph = graph

    def connect_to(self, other: Hashable, info: Any = None):
        """
        Connects this node to another node in the graph associated with this molecule.

        Notice how the other value does not have to be an instance of :class:`Node`.

        :param other: The value of the node to connect to
        :type other: Hashable
        :param info: The value associated with the new edge, defaults to None
        :type info: Any
        """
        self.graph.connect(self, other, info)

    def is_connected_to(self, other: Hashable) -> bool:
        """
        Determines if this node is connected to another node.

        Like :meth:`Node.connect_to`, the other node value
        does not have to be an instance of :class:`Node`.

        :param other: The value of the node to check
        :type other: Hashable
        """
        self.graph.are_connected(self, other)

    def neighbours(self) -> Dict[Hashable, Any]:
        """
        Returns a dictionary of the neighbours of this node,
        mapping the node to the info associated with the edge.

        :return: The dictionary
        :rtype: Dict[Hashable, Any]
        """
        return self.graph.neighbours_of(self)
