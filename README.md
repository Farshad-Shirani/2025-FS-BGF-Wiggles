# 2024-AmNat-FS-BGF

# Stability of Species’ Range Limits Set by Interspecific Competition

This repository contains MATLAB codes used in the article:

"Environmental “Wiggles” as Stabilizers of Species Range Limits Set by Interspecific Competition"

All codes are written in MATLAB R2021a.

## Description
These MATLAB codes numerically compute solutions of systems of partial differential equations (PDEs) that simulate adaptive evolution of geographic ranges of biological species in a one-dimensional geographic space.

## How to run the simulations:

- The code associated with each of the simulations in the abovementioned paper are provided in separate folders. The core of the code in different folders is the same, only some changes are applied to the `initializeSimulation` and `Main_Simulation` files, based on the specifics of each simulation. In particular, all codes included in the subfolder named `AuxiliaryFunctions` are the same for all simulations. These codes are associated with the numerical scheme that is used to solve the PDEs, and should not be changed, unless absolutely needed.

- Each folder contains a function named `initializeSimulation`. Model parameters, simulation parameters, and discretization parameters of the numerical scheme are first set by the user through this function. The model parameters are defined as global parameters. Therefore, the values of model parameters should only be changed through this function.

- Each folder contains a script named `Main_Simulation`, which is the main code that performs the simulations based on the parameters set in the function `initializeSimulation`. When a simulation is complete, the resulting solutions can be saved using the `save` command provided at the end of the script. The computed solutions are stored in a structure array named `populations`. The parameters are also stored in three different structures, named `modelParameters`, `simulationParameters`, and `discretizationParamaters`. The path provided to the `save` command should include the suborder name `Results`, so that all simulation results are organized in this subfolder. The `Results` folder is currently loaded with the simulation results that the authors have performed and used in the abovementioned paper.

- The results can be plotted using the script `Plotting` available in each folder. The path of the results to be plotted should be given in the `load` command at the beginning of the script. The parameter `incrementSize` inside the script should be set to the desired incremental time steps that the resulting curves should be plotted.

- The folders associated with climate change simulations also include an additional subfolder named `Data`. This subfolder includes the results obtained from another simulation that should be used in the current simulation for initializing the populations. 

## Which simulation (folder) is associate with which figure?

- The code in the folder `Limit Formation_Same` Species is associated with Figure 1 in the paper.

- The code in the folders `Competitive Exclusion` and `Marginal Coexistence` is associated with Figure 2 in the paper.

- The code in the folder `Stable Limit Formation_Upstream` is associated with Figure 3 in the paper.

- The code in the folder `Climate Change_Upstream` is associated with Figures 4 and 5 in the paper.

- The code in the folder `Climate Change_Downstreem` is associated with Figure 6 in the paper.
