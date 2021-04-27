# Twisted-loop
Monte-carlo simulator of shapes of twisted circular filaments in Timoshenko beam approximation.

# Package contents
This code can:
1. Generate energy-minimizing configurations for a given geometry and linking number.
2. Export shapes into bead format compatible with rigid molecule simulator:
 https://www.fuw.edu.pl/~rwaszkiewicz/rigidmolecule/

# References
- *Theory of Supercoiled Elastic Rings with Self-Contact and Its Application to DNA Plasmids* Bernard Coleman and David Swigon; J. of Elast. (2000)
- *Dynamics of an electrostatically charged rod in fluid* S. Lim, Y. Kim and D. Swigon; Proc. R. Soc. A, (2010).
- *Computation of writhe in modeling of supercoiled DNA* K. Klenin and J. Langowski; Biopolymers (2000).

# Contact
If you want to colaborate on some use of this code try starting here:

https://www.fuw.edu.pl/~rwaszkiewicz

# To compile
`g++ main.cpp -o main.exe`

# License
This project is licensed under MIT license

# Dependencies (included in repo)
- JSON parser is licensed under MIT license
- Eigen library is licensed under MPL2 license (used for matrix inversion)
- Args parser from stackoverflow by 0x90 and iain
