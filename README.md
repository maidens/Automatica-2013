# Automatica-2013
Supplementary materials for the paper: "Lagrangian methods for approximating the viability kernel in high-dimensional systems" by John Maidens, Shahab Kaynama, Ian M. Mitchell,  Meeko M. K. Oishi, and Guy A. Dumont



You will need the following MATLAB toolboxes installed and 
included on your path:

* The Multi-Parametric Toolbox (MPT)
* The Ellipsoidal Toolbox (comes with MPT)
* CVX


The following files are included:

* computeViab_Ellipsoid.m   - Implementation of the ellipsoidal algorithm
                            (Algorithm 2 in HSCC 2012 paper)
* computeViab_Polytope.m    - Implementation of the polytope method described
                            in Algorithm 3
* computeViab_SupportVect.m - Implementation of the support vector method
                            described in Algorithm 4
* run_time.mat              - Contains the data for Figure 7 (included 
                            because Figure_7.m takes a long time to run)
* grid_viab.mat             - Containts data for Figures 5 and 6 because
                            grid-based computation of viability kernel
                            takes a long time to run
* Figure_*.m                - Produces Figure * (for * = 4:9)
