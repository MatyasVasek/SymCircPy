# SymCircPy
[![pypi version](https://img.shields.io/pypi/v/symcirc.svg)](https://pypi.python.org/pypi/symcirc)
[![Downloads](https://static.pepy.tech/badge/symcirc)](https://pepy.tech/project/symcirc)
[![Powered by SymPy](https://img.shields.io/badge/powered%20by-SymPy-Green.svg?style=flat&colorA=gray&colorB=green)](https://www.sympy.org/)
![Tests](https://github.com/MatyasVasek/SymCirc/actions/workflows/tests.yml/badge.svg)
[![codecov](https://codecov.io/gh/MatyasVasek/SymCirc/branch/main/graph/badge.svg)](  https://codecov.io/gh/MatyasVasek/SymCirc)

A lightweight Python package for symbolic circuit analysis.

## Documentation

Documentation is available on [GitHub Pages](https://matyasvasek.github.io/SymCirc/).

---

## What can SymCircPy do?

### Symbolic and Semisymbolic Analysis
**SymCircPy** currently offers symbolic and semisymbolic **DC**, **AC** and **transient** small-signal circuit analysis.
It supports the following ideal circuit elements: **resistors, inductors, capacitors, independent sources, controlled sources, ideal operational amplifiers and coupled inductors with initial conditions**.
Transient simulation allows for initial conditions of capacitors and standard/coupled inductors.

### Spice-like netlist parsing
A custom netlist parser allows the user to input circuits in the standard SPICE format. Note that this has some quirks (e.g. only circuit elements, subcircuits and models are alowed). Nested subcircuits are not supported. The netlist format is modeled to work well with the GEEC online circuit simulator.

### Linearized transistors
**BJT, MOS and diode** can be used in AC and TF analysis. The operating point has to be provided externaly. The model is constructed using the parasitic values provided.

### Circuit model approximation 
> **Work in progress**

After parsing a netlist into a `Circuit` object, the `Circuit.approximate()` method can be used to approximate a large circuit into a simpler one according to the error control parameters.

### Multiple construction methods and solvers
 - **Tableau construction method**: to solve for all the voltages and currents in the circuit. (slower)

 - **Two Graph Modified Nodal**: solves for nodal voltages and select currents that are not easy to derive after solving. Every other result can be obtained using methods of the `Analysis` class. (faster)

 - **Sympy solve**: fast solve with the final results in a simplified form. Fully symbolic analysis is generally limited to a little over 10 symbolic elements due to symbolic explosion of the results. This is a well known limitation of symbolic analysis.

> **Work in progress**
 - **LED-DDD solve**: a custom layered expansion determinant decision diagram solver based on the work of **Shi, Guoyong, et al.** Very fast and capable of processing much larger circuits. This is an experimental feature included mostly for future experimentation with circuit approximation.




---

## Install
Use this command to install via pip
```
pip install symcirc
```
## Hard dependencies
### SymPy
**SymCirc** is a light-weight package. It needs only the [SymPy](https://www.sympy.org/) package, which is used for computations. It should get automatically installed with **SymCirc** when installed using pip. 
## Optional dependencies
### gmpy2
To achieve better performance install the [gmpy2](https://gmpy2.readthedocs.io/en/latest/) package. If **gmpy2** is installed, **SymPy** automatically uses it for integer operations. It significantly impacts semi-symbolic analysis performance. 
### symengine
[SymEngine](https://symengine.org) is an optional symbolic core for **SymPy**. To use the **SymEngine** core, run your script with the environment variable ```USE_SYMENGINE=1```, or run **AnalyseCircuit** with an optional argument ```use_symengine=True```
Example: ```symcirc.AnalyseCircuit(netlist, use_symengine=True)```
### numpy and matplotlib
For bode plots and other graphing utilities.
## Build

Or you can build with this command:
```
python setup.py sdist bdist_wheel
```

---

## Examples
See [examples](examples) for more insight into circuit analysis with SymCirc.
## Basic simulation example

``` python
from symcirc import *

# Insert your netlist
netlist = """CIRCUIT NAME - First line is always the circuit name
* This is a comment
R1 1 0 2k
R2 3 0 (1/G)
V1 2 1 dc 1 ac 1
R3 3 4 6k
C1 3 4 1n
R4 4 0 10k
V2 4 0 dc 5
I1 3 2 dc 1m ac 0
"""

# Execute netlist simulation
analysis_type = "DC"  # or "AC", "TF", "tran"
method = "tableau"  # Default is "two_graph_node", which is faster, but currently lacks coupled inductors.
symbolic = True  # If set to False, only elements which have no numeric value are left as symbolic. In this case only R2 stays symbolic.

circuit = AnalyseCircuit(netlist, analysis_type, symbolic=True, method=method)

all_values = circuit.component_values("all")
latex_formatted_values = to_latex(all_values)

print(all_values)
```

## Netlist syntax
There are many example netlists available in the [netlists](tests/netlists) folder. This folder contains a database of netlists which are used to automatically test every new version of SymCirc.

### Elements
* Resistor:          **RXXX N+ N-** _RESISTANCE_
* Capacitor:         **CXXX N+ N-** _CAPACITY IC=<INIT_VOLTAGE>_
* Inductor:          **LXXX N+ N-** _INDUCTANCE IC=<INIT_CURRENT>_
* Coupling:          **KXXX LYYY LZZZ K**
* Voltage source:    **VXXX N+ N-** _dc VOLTAGE ac AMPLITUDE <PHASE>_
* Current source:    **IXXX N+ N-** _dc CURRENT ac AMPLITUDE <PHASE>_
* IOAMP:             **AXXX Nout1 Nout2 Nin+ Nin-**
* VCVS:              **EXXX N1 N2 Ncontrol+ Ncontrol-** _GAIN_
* CCCS:              **FXXX N1 N2 VSENSE** _GAIN_
* VCCS:              **GXXX N1 N2 Ncontrol+ Ncontrol-** _GAIN_
* CCVS:              **HXXX N1 N2 VSENSE** _GAIN_

### Units
* terra: T
* giga: G
* mega: meg
* kilo: k
* mili: m
* micro: u
* nano: n
* pico: p
* femto: f




