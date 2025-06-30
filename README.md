# TMatrixCorrelationTools

## Description
This is a simple collection of MATLAB functions for correlation analysis using the transformation matrix **T**.  
No installation is required — just add the folder to your MATLAB path.

## Main Functions
- `T_Mass.m`   : Computes the T-Mass indicator to assess correlation in terms of mass between two models.  
- `T_Stiffness.m`: Computes the T-Stiffness indicator to assess correlation in terms of stiffness between two models.  
- `Rotmac.m`   : Computes the MAC and ROTMAC correlation matrices between two models.

## Usage Examples
Three example scripts are provided in the `Examples/` folder, demonstrating how to use the functions and interpret 
the results using a simple analytical 3-DOF system. For a deeper understanding, users are strongly encouraged 
to consult the related literature cited below.

## License
This code is released under the MIT License.  
You are free to use, modify, and distribute it, with proper attribution.

## Citing this work
If you use these functions in your research, please consider citing the following works:

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