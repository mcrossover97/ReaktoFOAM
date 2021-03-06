# ReaktoFOAM
A Geochemical Reactive Transport Solver for Acidizing: Coupling OpenFOAM and Reaktoro.

# How To Use
The ReaktoFOAM simulator consists of two OpenFOAM solvers: "LEAFOAM" and "PLEAFOAM" that solve the reactive transport problem using Local Equilibrium Approach and Partial Local Equilibrium Approach, respectively. To use the solvers, it is necessary to first install OpenFOAM:
https://openfoam.org/download/

And its plugin, swak4Foam:
https://openfoamwiki.net/index.php/Contrib/swak4Foam

Then, the Reaktoro library should be installed:
https://reaktoro.org/

After installing the Reaktoro library, it should be linked against the LEAFOAM/PLEAFOAM solver. For doing this, the user should include the 'include' and 'lib' directories and modify the 'option' files in the 'Make' directory of the solvers, accordingly:
https://www.cfd-online.com/Forums/openfoam-programming-development/231072-adding-libraries-openfoam.html

Please note that the PLEAFOAM solver considers the reaction rates. As a result, additional headers to describe the rates are found in the PLEAFOAM folder. ('ReactionRateFunctions.H' and 'reactions.H').

Should you have any question, please contact the developer: mcrossover97@gmail.com
