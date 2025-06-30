## Description
This folder contains example scripts to demonstrate how to use the correlation functions and how to interpret the results.

## Included scripts

- `Example1_StiffnessChange.m`  
  Simulates stiffness discrepancies with well-separated modes.   
  No example is included for stiffness discrepancies with closely-spaced modes, as similar results are expected.

- `Example2_MassChange_Separated.m`  
  Simulates mass discrepancies with well-separated modes.

- `Example3_MassChange_CloselySpaced.m`  
  Simulates  mass discrepancies with closely spaced modes.

## Theoretical summary

Users are strongly encouraged to consult the related literature before interpreting the results, as it provides the necessary 
theoretical background and practical context. However, a short summary is provided here for convenience.  

This code evaluates correlation between two modal models using a transformation matrix **T = pinv(B)·A**, where **B** and **A** are 
the modal matrices of systems B and A, respectively.  

Based on Structural Dynamic Modification (SDM) theory, the interpretation of the transformation matrix T depends on whether the observed discrepancies arise from changes in mass or stiffness.
Accordingly, the following can be stated:

|                        |&nbsp; Stiffness (Well-separated) &nbsp;|&nbsp; Stiffness (Repeated) &nbsp;|&nbsp; Mass (Well-separated) &nbsp;|&nbsp; Mass (Repeated) &nbsp;|
| :---                   | :---:                                  | :---:                            | :---:                             |:---:                        |
|T interpretation &nbsp;| `Tᵀ = R`                               | `Tᵀ = R`                         | `Tᵀ = R Tsh Tsc`                  | `Tᵀ = R Tsc `               |
|Rotation                | **✓**                                  | **✓**                           | **✓**                             | **✓**                      |
|Shear                   | ✘                                      | ✘                               | **✓**                             | ✘                          |
|Scaling                 | ✘                                      | ✘                               | **✓**                             | **✓**                      |

Based on the previous information, the following correlation indicators are derived:  
- `T_Mass`: Computes the angles between the column vectors of matrix **T**.  
- `T_Stiffness`: Computes the angles between the columns of matrix **T** and **Ω<sub>B</sub>²  T** or the angles between the columns of matrix **Ω<sub>B</sub> T**.  
- `Rotmac`: Computes the MAC between modal matrix **A** and the rotated modal matrix **B** (**B<sub>R</sub> = B · A**). It is a detector of shear effects.

## Related literature

1. **García Fernández, N.**, Fernández Fernández, P., Brincker, R., & Aenlle López, M. (2024).  
   *Mass and Stiffness Correlation Using a Transformation Matrix*, Infrastructures, 9(6), 96.  
   [https://doi.org/10.3390/infrastructures9060096](https://doi.org/10.3390/infrastructures9060096)

2. **Aenlle, M.**, García Fernández, N., & Fernández, P. (2024).  
   *Rotation of mode shapes in structural dynamics due to mass and stiffness perturbations*,  
   Mechanical Systems and Signal Processing, 212, 111269.  
   [https://doi.org/10.1016/j.ymssp.2024.111269](https://doi.org/10.1016/j.ymssp.2024.111269)

---

**This code was implemented by Natalia García Fernández.**  
For inquiries or suggestions, please contact: garciafnatalia@uniovi.es
