# EADYRKINST
MATLAB code for calculating the eigenmode solutions of the Eady model of Rossby-Kelvin instability formulated in Zurita-Gotor and Held (2025)

The dispersion relation and the various eigenmodes can be calculated by calling the function eigensolve as follows:

sweep=eigensolve(.01:.01:1.,-4:.1:8,-2:.05:2,@eadymodel,{3.5,-0.5});

This also requires the function eadymodel.m to be available. The function eigensolve is general and can be used to calculate the
eigenvalues in other configurations than that discussed in the paper, by simply calling eigensolve using the alternative model.

The data file control-lr.mat contains the output resulting from the above execution of eigensolve.

Type help eigensolve for additional details of its use
