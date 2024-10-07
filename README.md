# AUV_MPC

## Overview
This repository contains the code for the paper on Model Predictive Control (MPC) for Autonomous Underwater Vehicle (AUV) navigation. The code is designed to simulate and control the navigation of an AUV using MPC techniques.

## Installation
To get started with this project, clone the repository:

git clone https://github.com/Mohsen463/AUV_MPC.git
cd AUV_MPC

Ensure you have MATLAB installed on your system.

## Usage
To run the simulations and experiments, open MATLAB and navigate to the project directory. Then, run the main script:
1. Run the AUV_MPC_P2.m file, which is the main file. It calls the OceanEnvironment file and relevant functions to create the Ocean environment.
3. The system parameters can be changed in this file. It creates MPC matrixes and runs the P2 optimization problem (based on the paper) to initialize the paths.
4. Then it calls and runs P3, P4 and P5 optimization problems (described on the paper)1. Put all files in the same path or add paths.

## Contributing
Contributions are welcome! If you have any suggestions or improvements, please create an issue or submit a pull request.

## License
This project is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License. See the LICENSE file for details.

## Contact
For any questions or inquiries, please contact Mohsen (m.eskandari@unsw.edu.au).
Thank you for being so interested in our paper and the codes. 


