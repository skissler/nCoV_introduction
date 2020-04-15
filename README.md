# Simulation code associated with "Projecting the transmission dynamics of SARS-CoV-2 through the post-pandemic period"
S.M. Kissler*, C. Tedijanto*, E. Goldstein, Y.H. Grad+, and M. Lipsitch+

`*` denotes equal contribution

`+` denotes co-senior author

doi: 10.1126/science.abb5793

Correspondence: skissler@hsph.harvard.edu

__figuremaker.R__ reproduces Figure 2 in the manuscript and demonstrates how to import and manipulate the NREVSS data.

__invasion_rates.R__ is used within figuremaker.R for the ode simulation

__license.txt__ contains the GNU Version 3 license for use and redistribution of this code.

__nrevssCDC_ILI.csv__ contains the NREVSS data on the human betacoronaviruses, used in figuremaker.R

__postpandemicscenarios_full.nb__ is a _Mathematica_ file containing the simulation code for Figures 3-6 in the manuscript. It contains code for simulating the three-strain simulation model and the SARS-CoV-2 model with hospitalization and critical care arms. 

All code is available under the GNU General Public License, version 3, included in this repository under the following terms: 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
