# CPU optmization of Time Step Kernels in DEM Simulations 

- This piece of code is extracted from the [LIGGGHTS library] (https://github.com/CFDEMproject/LIGGGHTS-PUBLIC/blob/master/src/fix_check_timestep_gran.cpp)
which is a Discrete Element Method package used for simulations of granular flows.
- The class method fun() in this code implements simplified versions of two specific approaches to calculate estimates of time steps 
required for explicit time integration employed in the package.
- The performance (runtime for computation of the time step estimates) of the original unoptimized code is measured for various problem sizes 
(signifying the number of particles in the simulation). 
- We then perform two level of optimizations on this code: 
   - Level 1 optimizations corresponds to simply using GCC optimization flags (-O3 -march=native)
   - Level 2 optimizations corresponds to code refactorings combined with proper usage of GCC optimization flags (-O3 -march=native)
- Table shows performance gains at each optimization level (***About 13X speedups***)
 <div align="center">
   <img width="500" alt="cpuOptimizationTable" src="https://github.com/user-attachments/assets/59df7922-1ab1-4d7e-bfc3-035166697670"> 
 </div>
  
- Bar graph shows the same performance gains at each optimization level more visually (the linear scale-up in performance in visible)
<div align="center">
  <img width="423" alt="comparisonHistogram" src="https://github.com/user-attachments/assets/86b342a5-cee7-4d5a-8517-fa5ad78be9db">
</div>
