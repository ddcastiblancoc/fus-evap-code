# fus-evap-code
This Python code calculates some physical quantities and all possible residues related to a fusion-evaporation nuclear reaction with semi-classical calculations. This code is part of a B.Sc. Physics dissertation of the author and published in the Momento journal of the Universidad Nacional de Colombia as .

## Details
This article presents a study of fusion-evaporation nuclear reactions. Starting from a detailed description of the semi-classical theoretical framework behind this nuclear reaction, quantities such as the cross section of compound nucleus formation and various evaporation residues after its formation, as well as their cross sections (proportional to the events number), were estimated by means of a Python code. The code splits the compound nucleus formation process and its subsequent decay into several residual nuclei, which occurs as a sequential particle emission. In order to prioritize a first approximation theory, different nuclear models, with semi-classical and statistical origin, related to projectile-target fusion, light particle evaporation (n, p, $\upalpha$) and fission, were described in detail.

The purpose of this Python code is calculate cross sections of each exit channel of a fusion-evaporation reaction, given a projectile and a target.
This code uses Bass model and Wong formula (or classical) to calculate fusion cross section, Monte Carlo to select a decay mode randomly, Weisskopf-Ewing formula to calculate evaporation decay widths and Bohr-Wheeler model to calculate fission decay width.

## Using the code
It receives as input: mass number of projectile, atomic number of projectile, mass number of target, atomic number of target, laboratory energy of projectile, and number of cascades.
Compiling this code returns (in console): (1) data about projectile, target and compound nucleus, (2) physical quantities related to the formation of compound nucleus, (3) list of found residues from compound evaporation (and each evaporation cross section).

Please feel free to send me any question or suggestion to my email [ddcastiblancoc@unal.edu.co](mailto:ddcastiblancoc@unal.edu.co).

## References
