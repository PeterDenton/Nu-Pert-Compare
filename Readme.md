### Compiling the main code
1. Install libeigen3-dev (necessary for `Diagonalization.h`)

   You may also have to upgrade to the latest version via a ppa: https://launchpad.net/~nschloe/+archive/ubuntu/eigen-backports

2. Comment out parts of `main.cpp` as desired.
3. Run `make` in the main directory. `make clean` cleans the compilation files.
4. Run `./main`.
5. Numerical outputs are stored in `data/`. Python scripts for generating figures can be found in `py/`. The resultant figures appear in `fig/`.

This code has been tested on ubuntu and should work on other debian systems.

### Vacuum parameters
1. The vacuum oscillation parameters are given in `Parameters.h` and `Parameters.cpp`.
2. After modifying one of: `s12sq`, `s13sq`, `s23sq`, `Dmsq21`, `Dmsq31`, or `delta`, run `Recalc_Parameters();` to set `c12sq`, etc.
3. Mass squared differences are in eV<sup>2</sup> as usual and the raw angles/phases given by `t12`, `t13`, `t23`, and `delta` are all stored in radians.

### Function parameters
1. Each `Pmue` function is a function of a, L, and E, in that order (the order is easy to remember because of ale).
Some functions (DMP and AM) also carry an order parameter, see the header files.
2. The matter potential a is defined as 2EV<sub>CC</sub> and carries units of eV^2.
For example, given Y<sub>e</sub>=0.5 and &rho;=3.0 g/cc, one would set a by writing `a = E * 0.5 * 3.0 * YerhoE2a;` where YerhoE2a corrects the units.
3. The baseline and energy are given in km and GeV.
4. To reverse to &nu;<sub>e</sub> to &nu;<sub>mu</sub> oscillations take L to be negative.
5. For anti-neutrinos take E to be negative.
6. The choice of notation of the various expressions roughly follows those from the original papers which may differ from the associated article.

### Incorporating in your own code
1. The code for each expression is designed to be simple and stand alone.
2. Put the relevant cpp and h files in the appropriate locations for a given expression

   libeigen3-dev (see above) is necessary only for diagonalization

3. Include the parameters files (see above) and run `Recalc_Parameters();` after changing the vacuum parameters.

### Adding a new expression
1. Copy the template from, say, `Zeroth.h` and `Zeroth.cpp`. Modify the header guard and namespace. Insert your formula into the new cpp file.
2. Add an enumeration to the list in `Probabilities.h` anywhere before `LAST`.
3. Fill in the `Pmue` and `Name` functions with the new enum and include the header file in `Probabilities.cpp`.

### References
1. This code is free to use and modify.
2. If you use any aspect of this code or any related results, please reference the associated article **[arXiv:xxxx.xxxxx](https://arxiv.org)** as well as this code itself on zenodo **[doi:xxxx](https://zenodo.com)**.
3. If you use any of the included expressions, please also reference the article from which the expression comes from, each of which is mentioned in the relevant header file.

### Disclaimer
Every effort has been made to make the implementation of each expression as fast as possible while remaining true to the original expression.
If you find room for improvement in the implementation of an expression, please submit a pull request and it will be reviewed or add an issue.

