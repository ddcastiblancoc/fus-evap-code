# fus-evap-code
This Python code calculates some physical quantities and all possible residues related to a fusion-evaporation nuclear reaction with semi-classical calculations. This code is part of an undergraduate thesis of the author at Universidad Nacional de Colombia and published in an article form in the MOMENTO physics journal of the Universidad Nacional de Colombia[^1].

## Details
The structure of this code is based in a study of fusion-evaporation nuclear reactions, from a semi-classical perspective. Quantities such as the cross section of compound nucleus formation and various evaporation residues after its formation, as well as their cross sections (proportional to the events number), are the main objective of the calculations. The code splits the compound nucleus formation process and its subsequent decay into several residual nuclei, which occurs as a sequential particle emission. In order to prioritize a first approximation theory (much simpler to understand for new students interested in this topic), different nuclear models, with semi-classical and statistical origin, related to projectile-target fusion, light particle evaporation (n, p, $\upalpha$) and fission, were implemented here.

Then, given a target and projectile nucleus (at a certain energy), this code is capable of is calculate cross sections of each exit channel of this type of reaction. The theoretical framework used focuses in the Bass model and Wong formula (classical formula is also available) to calculate fusion cross section, Monte Carlo algorithm to select a random decay mode, Weisskopf-Ewing formalism to calculate evaporation decay widths and Bohr-Wheeler model to calculate fission decay width.

## Contents
This repository contains:
- A copy of the thesis on which is based this code (in spanish), in the file `Thesis (in Spanish).pdf`.
- The list of all routines used in the main code, in the file `routines_fus_evap_code.py`.
- The main code `main_fus_evap_code.py` which makes the calculations and shows all results in console.
- A list of the information of all isotopes found in NuDat database[^2] (extracted on Jan-22-2022) in `isotope_data_NuDat.txt`.

## Using the code
To run the code, you need all `.py` files to be in a single folder (and that has to be the working directory). Then, the main code can be compiled normally in your compiler of preference. The code uses common Python libraries like `numpy`, `pandas`, `scipy`, `random`, `joblib`, `multiprocessing`, `collections` and `time`, so make sure those are installed. It needs an additional library which needs to be installed manually via `pip`: `fisbar` by G. Henning[^3].

The code needs the input information of a reaction of interest, and outputs the results of the reaction in console. Especifically, it receives as input: mass number of projectile, atomic number of projectile, mass number of target, atomic number of target, laboratory energy of projectile, and number of cascades. It returns: (1) data about projectile, target and compound nucleus, (2) physical quantities related to the formation of compound nucleus, (3) list of found residues from compound evaporation (and each evaporation cross section).

Please feel free to send me any question or suggestion to my email [ddcastiblancoc@unal.edu.co](mailto:ddcastiblancoc@unal.edu.co).

[^1]: Castiblanco, D. (2024). Cross section estimation for heavy ion nuclear reactions with a cascade code of fusion evaporation. MOMENTO, (68), 52–68. [https://doi.org/10.15446/mo.n68.106810](https://revistas.unal.edu.co/index.php/momento/article/view/106810)
[^2]: https://www.nndc.bnl.gov/nudat3/indx_sigma.jsp.
[^3]: https://gitlab.in2p3.fr/gregoire.henning/fisbar-python/-/tree/v001.
