# AUV_MPC

Thank you for being so interested in our paper and the codes. 

1. Put all files in the same path or add paths.
2. Run the AUV_MPC_P2.m file, which is the main file. It calls the OceanEnvironment file and relevant functions to create the Ocean environment.
3. The system parameters can be changed in this file. It creates MPC matrixes and runs the P2 optimization problem (based on the paper) to initialize the paths.
4. Then it calls and runs P3, P4 and P5 optimization problems (described on the paper)
