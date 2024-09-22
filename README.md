## R2F2: a *R*un-time *R*econfigurable *F*lexible *F*loating Point Multiplier

### Paper: Exploring and Exploiting Runtime Reconfigurable Floating Point Precision in Scientific Computing: a Case Study for Solving PDEs

#### Cong Hao (Georgia Tech), ASP-DAC 2025

### Files:
- templated_floating_point.cpp
    - Arbitrary precision floating point multiplier, whose decision is decided at compile time
    - Acronym in code: **TFP**
- flexible_floating_point.cpp
    - The proposed R2F2, run-time reconfigurable flexible floating point multiplier
    - Acronym in code: **FFP**
- tb.cpp
    - The testbench for validating the accuracy of TFP and FFP
- heat_equation.cpp + plot_heat_equation.ipynb
    - C++ implementation of 1D heat equation and its visualization
- shallow_water_equation_python.ipynb
    - Original python implementation of shallow water equation (SWE)
- shallow_water_equation.cpp + plot_shallow_water_equation.ipynb
    - C++ implementation of SWE and its visualization
- script.tcl
    - Synthesis the designs into FPGA implementation