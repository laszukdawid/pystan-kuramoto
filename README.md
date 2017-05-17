# Kuramoto solver

Project intends to fit Kuramoto to provided data. 
Currently there are two available scenarios:
* **all-to-all**, where fit is for _N_ separated oscillators, and
* **all-to-one**, where fit is for signal that implicitly is composed of _N_ oscillators.

## General use
In all cases the model is ran from respective directory with
```
$ python kuramoto.py
```
If this is run for the first time it might take a while to execute as Stan needs to compile the model.
Another script is provided for displaying and comparing reconstruction:
```
$ python plot.py
```

Parameters for the models are stored in `model.config`. Currently this is used for generating model for tests and as initial values.

## Fit all-to-all
Input is supposed to be of dimension <i>N</i>x<i>T</i>, i.e. _N_ oscillators and their _T_ phases in time or Y<sub>N</sub>(t) = (y<sub>1</sub>(t), y<sub>2</sub>(t), ..., y<sub>N</sub>(t)).
The model fits initial positions, couplings and intrinsic frequencies.


## Fit all-to-one
The input is a linear superposition of many oscillators: S<sub>N</sub>(T) = s<sub>1</sub>(t)+s<sub>2</sub>(t)+...+s<sub>N</sub>(t).
The model fits initial positions, couplings and intrinsic frequencies.
