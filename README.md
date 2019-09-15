# Schubert-Polynomial-Package
A Mathematica package for studying the combinatorics of Schubert polynomials, permutations, and more. To run this package from another Mathematica notebook, place it any of the locations in your computer stored in the Mathematica path variable $Path and run <<SchubertPolynomials.m in your notebook. 

Alternatively, place it in "full_file_location" and add the following code to your notebook:

Block[{packages},
  packages = {"full_file_location"};
  If[! SubsetQ[$Path, packages], $Path = 
    DeleteDuplicates[Join[$Path, packages]]];
  << SchubertPolynomials.m;
  ];
