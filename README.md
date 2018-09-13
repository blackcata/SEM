# SEM (Synthetic Eddy Method)
Inflow generation method by using Synthetic Eddy Method(SEM) which is extended to temperature not only fluid components.
At first, this method is developed by Jarrin (Jarin et al.2005) with a basic idea that turbulence is a superposition of coherent structures, 
which are called eddies. This code is extended to include not only fluid components but temeperature too. 
Code is modulized to apply any model easily and parallelized by OpenMP to make acceleration of compuational time.
Also various type of fluids are tested to validate this code's performance and their results are attached on related paper below.

This code's algorithm is applied to boundary layer and now sumiited to International Journal of Heat and Mass Transfer.

#Title: Extended Synthetic Eddy Method to Generate Inflow Data for Turbulent Thermal Boundary Layer
Listed co-author(s): Geunwoo Oh Eugene KyungMin Noh,
Corresponding Author: Jungil Choi

### 1. Tested fluid 
Each tested datas are saved at DATA folder. 
If you want to test or apply another type of fluid, 
you must make input file as prescribed slice file on DATA folder - 2D slice data including mean terms and Reynolds stress terms including temeperatue components. 
  - Channel case (Re_tau = 180, Ri = 0) - Garcia et al.2011 
  - Turbulent boundary layer case - Kong et al.2000
  - Stratified stable boundary layer case - Basu.2006
  
### 2. Related papers  
  - Related paper : https://goo.gl/eEgVpV
  - Related report : https://goo.gl/4B24AC

### 3. Preliminary setting for SEM
  - Set the number of meshes (Ny, Nz - SETUP subroutine)
  - Set the number of eddies and eddy length scale (N,SIGMA - this'll be modified to automatically calculate soon)
  - Set the number of iteration time (Nt - flow_module)
  - Set the number of OMP threads
  
### 4. Ouput files
  - Mean_profiles.plt
  - RMS_profiles.plt
  - Edd_POS.plt
  - U_ins.plt
  
### 5. Code composition 
Each codes consists of each subroutine and module and detail explanation of code and variables are described at each code. 
  - SEM_main.f90
  - SEM_module.f90
  - SEM_setup.f90
  - SEM_read.f90
  - SEM_eddy_setting.f90
  - SEM_fluctuation_gen.f90
  - SEM_combine_slice.f90
  - SEM_convection.f90
  - SEM_statistics.f90
  - SEM_write.f90
  
  At module_seperation & SEM_SBL branch, each codes are modulized to easily combine to other flow solvers without any dependencies.
  - main.f90
  - SEM_module.f90
  - IO_module.f90
  - flow_module.f90
  - SEM_write.f90
    
'main.f90' and 'flow_module.f90' contents has to be replaced by flow solver properties and codes to 'SEM_module.f90' to applied this SEM.
After validation of the fluid data which want to test, remove SEM_STAT and SEM_write subroutine.
  
  
