## Data for the paper, "Transition Energy Analysis of C15-like Defect Clusters in Tungsten using Molecular Dynamics Simulations"

#Purpose:
This document describes all data, code, and LAMMPS input files required to 
reproduce the results of the manuscript.

#Folder Structure
    * Data/
        - combined_rampflip_data.json
        - combined_transflip_data.json

        
    * Code/
        - ramp.py
        - Isothermal.py

        
    * LAMMPS/
    
        *[NPTeqlb]/
        *Cluster type/
        *Name/
        -vis.dat  (represents the embedded C15-like defect structures for NPT simulations.)
        *[eqlb]/
        - eqlb.in (Input files used for equilibration runs of all eleven C15-like defect clusters) 
        
        *[Ramp]/
        *Cluster type/
        *Name/
        *slowramp/
        - ramp.in (Input file used for ramped-temperature simulations of all eleven C15-like defect clusters)
        
        *[Trans]/
        *Cluster type ConstT/
        *Clsuter type Name/
        *Name/
        *ConstT/
        *Temp/
        -Trans.in (Input file for isothermal transition simulations of the selected clusters: 2C, Tripod-3B, and 5C)

## Potentials Used

        - Due to copyright restrictions, the potential file cannot be redistributed as part of this repository. However, it is publicly available from the **OpenKIM** repository:
         https://openkim.org/id/EAM_MagneticCubic_DerletNguyenDudarev_2007_W__MO_195478838873_002

        -The stiffening of this potential is discussed in the following reference: https://www.sciencedirect.com/science/article/abs/pii/S0168583X09007575
        
## How to Run

## 1. LAMMPS Simulations
   Raw molecular dynamics data can be generated using LAMMPS input files located in the `LAMMPS/` directory.
          -Example command for ramped-temperature simulations: lmp_serial -in ramp.in and for isothermal transition simulations: lmp_serial -in Trans.in
          
## 2. Python Analysis Scripts
   Run the executable Python scripts:
          - python ramp.py
          - python Isothermal.py
          - These scripts read the JSON files in the Data/ directory, generate all analysis figures, and perform Arrhenius and Eyring fits.



