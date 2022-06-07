# fus-evap-code
This Python code calculates some physical quantities and all possible residues related to a fusion-evaporation nuclear reaction with semi-classical calculations. This code is part of a undergraduate thesis of the author at Universidad Nacional de Colombia.

## Details
This article presents a study of fusion-evaporation nuclear reactions. Starting from a detailed description of the semi-classical theoretical framework behind this nuclear reaction, quantities such as the cross section of compound nucleus formation and various evaporation residues after its formation, as well as their cross sections (proportional to the events number), were estimated by means of a Python code. The code splits the compound nucleus formation process and its subsequent decay into several residual nuclei, which occurs as a sequential particle emission. In order to prioritize a first approximation theory, different nuclear models, with semi-classical and statistical origin, related to projectile-target fusion, light particle evaporation (n, p, $\upalpha$) and fission, were described in detail.

The purpose of this Python code is calculate cross sections of each exit channel of a fusion-evaporation reaction, given a projectile and a target.
This code uses Bass model and Wong formula (or classical) to calculate fusion cross section, Monte Carlo to select a decay mode randomly, Weisskopf-Ewing formula to calculate evaporation decay widths and Bohr-Wheeler model to calculate fission decay width.

## Contents
This repository contains:
- A copy of the thesis on which is based this code (in spanish), in the file `Thesis (in spanish).pdf`.
- The list of all routines used in the main code, in the file `routines_fus_evap_code.py`.
- The main code `main_fus_evap_code.py` which makes the calculations and shows all results in console.
- A list of the information of all isotopes found in NuDat database[^1] (especially for mass defects) in `isotope_data_NuDat.py`.

## Using the code
To run the code, you need all `.py` files to be in a single folder (and working directory). Then, the main code can be compiled in your compiler of preference.
It receives as input: mass number of projectile, atomic number of projectile, mass number of target, atomic number of target, laboratory energy of projectile, and number of cascades.
Compiling this code returns (in console): (1) data about projectile, target and compound nucleus, (2) physical quantities related to the formation of compound nucleus, (3) list of found residues from compound evaporation (and each evaporation cross section).

fisbar needs to be installed manually; see: gitlab.in2p3.fr/gregoire.henning/fisbar-python/-/tree/v001

Please feel free to send me any question or suggestion to my email [ddcastiblancoc@unal.edu.co](mailto:ddcastiblancoc@unal.edu.co).

[^1]: www.nndc.bnl.gov/nudat3/indx_sigma.jsp
