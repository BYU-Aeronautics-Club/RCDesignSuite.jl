# RCDesignSuite

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://aeronautics.byu.edu/RCDesignSuite.jl/stable)
![](https://github.com/BYU-Aeronautics-Club/RCDesignSuite.jl/workflows/Run%20tests/badge.svg)


A suite of tools for conceptual and preliminary design of basic RC airplanes based on public FLOW Lab codes.


## Dependencies:

<!-- - [Atmosphere.jl](https://github.com/byuflowlab/Atmosphere.jl) : a standard atmosphere model -->
- [CCBlade.jl](https://github.com/byuflowlab/CCBlade.jl) : a blade element momentum code with guarenteed convergence, suitable for gradient-based optimization
- [Composites.jl](https://github.com/byuflowlab/Composites.jl) : a compilation of composite structure theory calculations
- [FLOWMath.jl](https://github.com/byuflowlab/FLOWMath.jl) : various mathematical tools, including Akima splines.
- [GXBeam.jl](https://github.com/byuflowlab/GXBeam.jl) : a pure Julia implementation of Geometrically Exact Beam Theory
- [MotorPower.jl](https://github.com/byuflowlab/MotorPower.jl) : a compilation of various electric motor related models
- [VLM.jl](https://github.com/byuflowlab/VLM.jl) : a VLM code similar to AVL by Mark Drela, but suitable for gradient-based optimization
- [Xfoil.jl](https://github.com/byuflowlab/Xfoil.jl) : a bare-bones julia wrapper for xfoil (from pyxlight)


## Additional Capabilities

- Conceptual Design Tools : suite of tools for conceptual design of RC aircraft
- APC Propeller Model : input diameter and pitch and output all the required rotor information to run CCBlade