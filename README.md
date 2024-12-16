# SatAndUnsatPorom
Article: Modelling transition between fully and partially saturated porous media via a new implicit formulation of interface evolution


In this repository, one can find the four different code used to generate the data presented in the article.

BeautifiedUFL is the folder gathering all the UFL files required to run the computations.
Mesh folder holds mesh file used in the computations.
Ref, Sink, rate7-8 and pH120 are four different folders serving as example of results you obtain when running the script "script.cpp" laying inside.
Ref is the folder corresponding to the 1D test case, Sink corresponds to the 2D test case, rate7-8 corresponds to the speed rate sensitivity analysis for 175 time steps and pH120 corresponds to the test using 120% of the standard characteristic pressure of the model.

"pH120"'s script is made to be run with one argument x corresponding to the x % of standard characteristic pressure one would employ.
"rate7-8"'s script is made to be run with two arguments x and y corresponding to imposing 200x/y time steps of liquid pressure decrease/increase. 

The code requires FEniCS legacy version available at this adress:
https://launchpad.net/~fenics-packages/+archive/ubuntu/fenics

