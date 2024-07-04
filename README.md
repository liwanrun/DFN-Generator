# Two-dimensional Rough Discrete Fracture Network (RDFN) generator

## Introduction

In order to realistically simulate the morphology of natural fractures embeded in the engineering rock masses, a Fourier transform approach is utilized to controllably reconstruct the fracture surface roughness. Combined with Monte-Carlo sampling technique, Rough Discrete Fracture Networks (RDFNs) can be statistically recomplished, which plays a fundamental role in computational rock mechanics and rock engineering.

We intended to incorparate the reconstructed RDFN into our in-house code series QFDEM/YFDEM to study the problem of hydraulic fracturing, by means of the mesh-based methods including Finite Element Method (FEM), Finite Difference (FDM) and Finite Volume Method (FVM). The current repository is belong to the preprocessing module of further numerical analysis,

## Requirements

The functionality is implementd in Python, along with common third party libraris, which consists of [numpy](https://numpy.org/), [matplotlib](https://matplotlib.org/), [shapely](https://shapely.readthedocs.io/en/stable/manual.html) and [gmsh](https://gmsh.info/). You can easily install these pacakages like this:

```shell
pip install numpy, matplotlib, shapely, gmsh
```

## Results
