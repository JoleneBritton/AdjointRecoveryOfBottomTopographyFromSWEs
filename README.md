# AdjointRecoveryOfBottomTopographyFromSWEs
The code base corresponding the the paper 'Recovery of a Time-Dependent Bottom Topography Function from the Shallow Water Equations via an Adjoint Approach' by Jolene Britton, Yat Tin Chow, Weitao Chen, and Yulong Xing.

### To run examples from the paper:
- **SWE_Run.m:** this will generate the figures that come in triple (best iteration vs true p, various iterations, error vs iteration number plots). See the comments to understand which set of parameters corresponds to each figure.
- **LCurveTest.m:** this will generate the L-Curve plot. The data must already be generated for this to work (see appropriate cases in SWE_Run to run prior to executing this code)
- **AccuracyTest.m:** this will generate the errors and orders of accuracy seen in Table 3

### Other files:
- **SWE_Main.m:** this is the driver file that runs the entire iterative algorithm
- **SWE_Data.m:** stores the properties of the data class and contains the function that populates the data class instance
- **SWE_Meas.m:** generates the measured data using the true p function
- **SWE_Forward.m:** driver code for the forward problem using the RK3 DG method which computes the state variables for a given iteration p^k
- **SWE_Inverse:** contains 
    - the function for computing the update to p^k based on Algorithm 3.1
    - the driver code for the RK3 DG method for the adjoint problem which computes the solutions to the adjoint and the components needed to calculate the gradient
- **SWE_Update:** contains the calculations for the RK3 methods - the integral, flux, and source term calculations
- **SWE_SetUp:** contains a series of simple functions including
    - projection of ICs into polynomial space
    - determining cell interface values
    - computing boundary conditions
    - computing the mass matrix
    - computing the Lax-Friedriches coefficient
    - converting variables to a more coarse mesh
    - more
- **SWE_Limiters:** the slope limiter code needed for when discontinuities arise
- **SWE_PostProcessing:** contains 
    - the plotting functions that get used during the iterative process
    - the function to calculate errors for the accuracy test
- **SWE_Results:** generates all the final plots at the end of all the iterations, the plots used in the paper
- **lgwt.m & lglnodes.m:** used for nodes and weights in numerical integrals
