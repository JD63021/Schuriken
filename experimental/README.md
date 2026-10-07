# Experimental polynomial Schur / FGMRES work

The P10 polynomial-Schur FGMRES path from the research campaign is deliberately
not part of Schuriken v0.1.  It showed strong acceleration on some pipe cases,
but was more sensitive to under-relaxation and pressure tolerances than the
stock-operator + adaptive-POD path.  The benchmark CSV retains representative
results so this backend can be revisited later without making it the default.
