# TEMFit

## TEMFIT Overview

The TEMFIT package is a biophysics-based fitness model designed to determine the predicted bacterial fitness upon point mutations on TEM-1 B-lactamase. The package calculates relative fitness of different TEM-1 mutants at different antibiotic concentrations using experimentally measured kinetic rate constants. The model determines the concentration of B-lactam in the periplasm, and the relative fitness of mutants to the wild-type. 

Three input files are required to run TEMFIT:
1.	Cell specifications with kinetic and cellular constants for the cell.
2.	Enzyme specifications with kinetic variables for B-lactamase.
3.	Dosage specifications with specific antibiotic concentrations.

Depending on the selected efflux mode, TEMFIT has two running modes. It can assume equilibrium between diffusion and hydrolysis (efflux = 0), or it can calculate relative fitness by considering the diffusion and hydrolysis apart (efflux = 1).

The package first checks if the “deSolve” package is installed in the user’s R library, and loads it accordingly. TEMFIT then loads the variables and values specified in the three input files (fileReader.R). If the selected mode is efflux 0, then relative fitness is calculated directly. If the selected mode is efflux 1, the relative fitness is calculated by solving the differentials in TimeTEM.R. A system of equations is used to solve for the concentration of B-lactam in the periplasm, the concentration of B-lactam in the media, the concentration of peptidoglycan, and the concentration of TEM-1. These are solved by determining the diffusion rate, the peptidoglycan-binding-protein (PBP) inhibition rate, the hydrolysis rate, and the efflux rate. 
Finally, TEMFIT outputs the results in a file which contains the rates of change used to solve the equations, the concentration of B-lactam in the periplasm, and the relative fitness of mutants to the wild-type.

#Add link to article

#License
Copyright (C) 2018 Serohijos Lab

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
