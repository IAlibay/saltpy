saltpy
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/IAlibay/saltpy.svg?branch=main)](https://travis-ci.com/IAlibay/saltpy)
[![codecov](https://codecov.io/gh/IAlibay/saltpy/branch/main/graph/badge.svg)](https://codecov.io/gh/IAlibay/saltpy/branch/main)


A python toolset for calculating and adding salt concentration to solvated atomistic systems

## Contents

### estimators:
Several methods to estimate the number of ions to add to a solvated solute box.

These include:
- `add_neutralize` (the add-then-neutralize method)
- `genion` (gmx genion's version of add-then-neutralize)
- `sltcap` (as described by Schmit et al., JCTC, 2018)
- `split` (as described by Machado and Pantano, JCTC, 2020)

## calculate
An MDAnalysis driven method for using the estimators on an input solvated structure.

## ionize
Wrappers for common MD engines to add ions.

Currently implemented:
- `genion` (using GromacsWrapper)

### Copyright

Copyright (c) 2020, Irfan Alibay


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.
