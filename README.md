# Stochastic EM Channel Modeling
The stochastic EM channel modeling in scattering environment using a probabilistic framework. For the technical details please refer to the paper "**Electromagnetic Channel Modeling and Capacity Analysis for HMIMO Communications**" (has been accepted for publication in IEEE Transactions on Wireless Communications). 

>This paper presented a universal EM-compliant channel model for HMIMO communication systems based on stochastic Green’s functions. By incorporating a large number of planar waves, the proposed channel can be statistically represented
in a probabilistic framework. The asymptotic channel capacity and capacity region are derived based on the physics-oriented channel models.

<ins>One should cite this paper for using these codes.</ins>

MATLAB code prepared by Li Wei (e-mail: l_wei@ntu.edu.sg and weili_xd@163.com)

> [!NOTE]
> If you want to incorporate the impact of antennas (e.g., S-parameter, impedance matrix, or radiated efficiency) in the modeled channel, i.e., coupled EM channel modeling, the **CST** EM simulation software is required, we provide some data samples (e.g., 2*2 antenna array with spacing 0.4\lambda at 5GHz) in this file. For more antenna setting, please run **main_CST_CouplingMatrixGen.m** to generate the required data for further channel modeling (please ensure that you have CST software before running  **main_CST_CouplingMatrixGen.m**).

> [!TIP]
> You can run the code directly by commenting out the mutual coupling matrix part, i.e., uncoupled EM channel modeling, if you don't have CST software.
 

