#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="A program for parsing SMILES chemical formulae to skeletal formulae"
    )
    parser.add_argument(
        "formula", help="The formula to parse",
    )

    return parser, parser.parse_args()


def main():
    parser, args = parse_args()


if __name__ == "__main__":
    main()
