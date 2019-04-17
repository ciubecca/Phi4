## Intro

This repository contains several branches. The "main" is an old branch which contains code for the 2d problem. 
The latest code for the 3d problem is called "nlo"

## Compile

Some of the most critical parts of the code are written in Cython. In order to compile the Cython file, run the following command:

```bash
python setup.py build_ext --inplace
```

## Test

In order to validate changes made to the code, I created a test module that compares the results to expected output. It can be run with 
```bash
pytest testBasis.py
```

## Run

Schematically, the code is split into "core" files, "generation" files (to compute eigenvalues and store them in a sqlite database) and "plotting" files which retrieve eigenvalues from the database and plot them. 
First, create the "data" folder to store the database:
```bash
mkdir data
```
For example, the following sequence of commands computes eigenvalues for (L=5, ETmax=24, g2=0), (L=6, ETmax=20, g2=0), (L=7, ETmax=18, g2=0), for various values of g4, 
```bash
python eigsET.py 5 24 0
python eigsET.py 6 20 0
python eigsET.py 7 18 0
```
Finally we can plot those eigenvalues as function of ET for g4=1.0, for instance, with
```bash
plotvsET_L.py 1.0 0
```
In general, both in the "generation" and "plotting" files there are various hardcoded parameters that should be changed by hand.
