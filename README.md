This code contains a number of methods to compute important quantities such as the self-force, self-inductance, and magnetic field for HSX coils and circular coils with circular cross-sections. (It should be easily adapted for any coil with circular cross-section described with a Fourier representation). All formulas referenced are found or can be determined from those presented in (Hurwitz, Landreman, and Antonsen; 2023). 

The code also contains a number of routines in "data_visualization.cpp" to export data to text files for plotting purposes. The file breakdown is as follows:

-"b_calcs.cpp" contains hifi and lofi methods for calculating the magnetic field
-"f_calcs.cpp" contains hifi and lofi methods for calculation the self-force
-"l_calcs.cpp" contains hifi and lofi methods for calculation the self-inductance

-"integration.cpp" contains helper code for numerical integration using GSL
-"point.cpp" describes our custom Frenet-Serret coordinate system
-"vectors.cpp" contains basic vector algebra methods
-"wire.cpp" contains information on specific coils and gives important coil-dependent quantities such as the position vector and unit vectors


