# usadeleq
solving the usadel-eq with a relaxation method from numerical recipes

This repository contains code supporting the bachelor thesis of Tim Lüders, University of Greifswald.
As described in the thesis, the used code consists of the files 


  nr3.h         (provided by numerical recipes, under copyright)
  function.h    (given here)
  Difeq.h       (given here)
  Solvde.h      (provided by numerical recipes, under copyright)
  Solving.h     (given here)
  main.cpp      (given here)
  
  with minor additions to nr3.h and Solvde.h discussed in the thesis. 
  
  With this code and the additional file one can solve the Usadel-transport-equation and from this obtain the retarded green's function, which yields the local density of states.
  With the density of states we can deduct important properties of the present System, such as if there are superconducting effects anywhere over the domain.
