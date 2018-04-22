# fluid

Free Fluid Flux-Barriers Rotor for Synchronous Reluctance Motor Drawing

[![DOI](https://zenodo.org/badge/123190369.svg)](https://zenodo.org/badge/latestdoi/123190369)

<img src="doc/gfx/3b.png" width="256">

`SynRM`, `SyRM`, `SynchRM`, `REL`, `SynchRel`,
`PMASynRM`, `PMASyRM`, `PMASynchRM`, `PMAREL`, `PMASynchRel`


## Goal

Provide a ready-to-use fully parametric drawing of the
Synchronous Reluctance Rotor with *fluid* flux-barriers.

The scope of this project is the computation of the flux-barriers points.
The drawing scripts are for demonstration purposes only.


## Requirements

Matlab or Octave or Python (NumPy + SciPy) to compute the points. 
The points calculation is general, so it could be implemented in any
language, but I chose Matlab/Octave because it is my standard 
interface with [FEMM](http://www.femm.info) software.

If you do not use FEMM, you can still use the calculation part and
make a porting for your CAD engine or FEA software.
If you do so, consider contributing to the project adding your 
interface scripting.

