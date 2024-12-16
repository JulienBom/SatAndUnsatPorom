#!/bin/sh


ffc -l dolfin Dual2D.ufl
wait
ffc -l dolfin Sum2D.ufl
wait
ffc -l dolfin InitialField2D.ufl
wait
ffc -l dolfin InitialNewton2D.ufl
wait
ffc -l dolfin Newton2D.ufl
wait
ffc -l dolfin SatProj2D.ufl
wait
ffc -l dolfin SurfInt2D.ufl
wait
ffc -l dolfin SeaProj2D.ufl
wait
ffc -l dolfin QuasiStat2D.ufl
wait
ffc -l dolfin VolDef.ufl
wait
ffc -l dolfin NormalFlux.ufl
wait
ffc -l dolfin SubdomainMaker.ufl

wait
ffc -l dolfin SubdomainMaker4.ufl
wait
ffc -l dolfin SubdomainMaker5.ufl
wait
ffc -l dolfin SubdomainMaker6.ufl
wait
ffc -l dolfin SubdomainMaker7.ufl

wait
ffc -l dolfin VolDef.ufl
