PyFormulate
===============

A library for parsing SMILES (simplified molecular-input line-entry system) strings, and displaying molecules.

SMILES
------

PyFormulate follows the `OpenSMILES specification`_ in its parser. 

Recently, IUPAC have adopted the OpenSMILES specification as the basis for a new standard called `SMILES+`_.
At the time of writing, the SMILES+ specification is exactly the same as the OpenSMILES standard, however as SMILES+ develops,
it may become the standard followed by PyFormulate.

To-do list
----------

- [x] Syntax

  - [x] Atoms

    - [x] Bracket Atoms
    - [x] Aliphatic Organic Atoms
    - [x] Aromatic Organic Atoms

  - [x] Chirality

    - [x] Tetrahedral
    - [x] Cis-trans
    - [x] Tetrahedral allene-like
    - [x] Square-planar
    - [x] Trigonal bipyramidal
    - [x] Octohedral

  - [x] Hydrogen count
  - [x] Charge
  - [x] Atom class
  - [x] Bonds

    - [x] Bond symbols
    - [x] Ring bonds

  - [x] Chains

    - [x] Branched atoms
      
      - [x] Branches

  - [x] Terminator

- [-] Semantics

  - [x] Elements

    - [x] Valid bracket elements
    - [x] Valid bracket aromatics
    - [x] Valid organic aliphatics
    - [x] Valid organic aromatics

  - [x] Organic subset

    - [x] Implicit hydrogen count

  - [x] Chirality

    - [x] Valid chiral class
    - [ ] Interpreting @ and @@ to the appropriate exact class

  - [x] Hydrogen count

    - [x] Illegal hydrogen count for hydrogen

  - [ ] Aromatics

    - [ ] Valid aromatic rings
    - [ ] Determining bonding between aromatic atoms
    - [ ] Adding hydrogens to organic aromatic atoms

.. _`OpenSMILES specification`: https://opensmiles.org/opensmiles.html
.. _`SMILES+`: https://github.com/IUPAC/IUPAC_SMILES_plus
