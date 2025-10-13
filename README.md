# PartonShower.jl

[![Build Status](https://github.com/Zykerin/PartonShower.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Zykerin/PartonShower.jl/actions/workflows/CI.yml?query=branch%3Amain)

<h2> Description </h2>


A simple angular-ordered final-state QCD parton shower written in Julia [[1]](#1)[[2]](#2). This shower is valid for any event that has colorless initial-state particles.



<h2> Installation</h2>

This package does not require any special steps to install. Just run Julia in the terminal and type

```julia
using Pkg
Pkg.add("PartonShower")
```

<h2>Usage</h2>

To run the shower, you need to have a Les Houches Event (LHE) file to use as an input. An example of the proper syntax using the test file `eejj_ECM206.lhe.gz` is below.

```julia
showeredEvents, energy = ShowerLHE("eejj_ECM206.lhe.gz")
```

This returns the showered events under `showeredEvents` and the center of mass energy under `energy`. The events are returned as a list of `Event` structures. 

Along with the `Event` structure there is the `Particle` and `Jet` structures which are layed out below.

- `Particle` - Contains the information for one particle

    - `id` - The PDG id of the particle (ex. 21 = gluon).
    - `status` - The status code of the particle. 
    - `t` - The starting scale for this particle.
    - `z` - The light-cone momentum fraction of this particle.
    - `m` - The mass of the particle.
    - `pT` - The transverse momentum of the particle.
    - `px` - The x-component of momentum.
    - `py` - The y-component of momentum.
    - `pz` - The z-component of momentum.
    - `E` - The energy of the particle.
    - `phi` - The azimuthal angle of the particle around the z-axis.
    - `qT` - The remaining components of the momentum according to the Sudakov basis [[1]](#2)[[3]](#3). 
    - `color` - The color of the particle in the form xxx.
    - `anticolor` - The anticolor of the particle also in the form xxx.
    - `alpha` - The $\alpha$ value in the Sudakov basis.
    - `beta` - The $\beta$ value in the Sudakov basis.
    - `virtuality` - The virtuality of the particle. 
    - `children` - A list containing the particle's children, if any. 

- `Jet` - A structure for a singular jet of an event which contains:

    - `AllParticles` - A list of all the particles in the jet.
    - `Progenitor` - The particle that the jet results from.

- `Event` - A structure for one event which contains:

    - `AllParticles` - A list of all the particles in the event.
    - `Jets` - A list of the jets in the event.


The showered events can be used for further analysis, or can be written out to another LHE file using the `WriteToLHE` function. Using an example filename `Showeredeejj_ECM206.lhe`, and the showered events in the previous example, the proper syntax is:

```julia
sigma = 1.2
sstdev = 0.6
WriteToLHE(showeredEvents, "Showeredeejj_ECM206.lhe", energy, sigma, sstdev)
```

The `sigma` value is the evaluated cross section of the events and the `sstdev` is the error of it.  

<h2> References </h2>

<a id ="1">[1] </a>
Bahr, M., and Others. ‘Herwig++ Physics and Manual’. Eur. Phys. J. C, vol. 58, 2008, pp. 639–707, https://doi.org/10.1140/epjc/s10052-008-0798-9. arXiv.

<a id ="2"> [2] </a> 
Papaefstathiou, Andreas. ‘$\texttt{Pyresias}$: How To Write a Toy Parton Shower’. arXiv [Hep-Ph], 6 2024, arxiv.org/abs/2406.03528. arXiv.


<a id="3">[3] </a> 
Gieseke, Stefan, et al. ‘New Formalism for QCD Parton Showers’. JHEP, vol. 12, 2003, p. 045, https://doi.org/10.1088/1126-6708/2003/12/045. arXiv.