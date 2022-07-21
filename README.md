# JuliaSONG: Julia Simulation Of Neutrino and Gamma-rays
A Julia version of FIRESONG

## Usage
Define the universe with `JuliaSONG.Cosmology`.
It should contain the $\Lambda$-CDM parameters and the function that describe the evolution.
The evolution function need not be normalized. There are 2 provided, `JuliaSONG.SFR` and 
`JuliaSONG.NoEvo`. The SFR is from Madau et al. The default cosmology uses `JuliaSONG.SFR`.

Generate the sources with `JuliaSONG.generate_source`.
E.g., for a source population with local density = $10^{-7} \mathrm{Mpc^{-3}}$, 
and the diffuse neutrino flux = $1.01\times10^{-8}\times\frac{E}{100~\mathrm{TeV}}^{-2.28} ~\mathrm{GeV/cm^2/s}$,

```JuliaSONG.generate_source(1e-7, 1.01e-8, 2.28)```.

The output will be the fluxes and the redshifts of the sources, in julia DataFrame format. 