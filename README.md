# PartonShower.jl

[![Build Status](https://github.com/Zykerin/PartonShower.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Zykerin/PartonShower.jl/actions/workflows/CI.yml?query=branch%3Amain)

<h2> Description </h2>


A simple colored final state parton shower written in Julia. This shower is valid for any event that has colorless initial state particles.



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

    - `id` - What type of particle it is, (ex. 21 = gluon)
    - `status` - The status of the particle
    - `t` - The starting scale for this particle
    - `z` - The light-cone momentum fraction of this particle
    - `m` - The mass of the particle
    - `pT` - The transverse momentum of the particle
    - `px` - The momentum in the x-axis
    - `py` - The momentum in the y-axis
    - `pz` - The momentum in the z-axis
    - `E` - The energy of the particle 
    - `phi` - The azimuthal angle of the particle
    - `qT` - The remaining compontents of the momentum according to the Sudakov basis
    - `color` - The color of the particle of the form 5xx
    - `anticolor` - The anticolor of the particle also of the form 5xx
    - `alpha` - The $\alpha$ value in the Sudakov basis 
    - `beta` - The $ \beta $ value in the Sudakov basis
    - `virtuality` - The virtuality of the particle
    - `children` - A list containing the particle's children, if any

- `Jet` - A structure for a singular jet of an event which contains

    - `AllParticles` - A list of all the particles in the jet
    - `Progenitor` - The particle that the jets results from 

- `Event` - A structure for one event which contains

    - `AllParticles` - A list of all the particles in the event
    - `Jets` - A list of the jets in the event

Using the returned showered events, you can either use them as you wish or write them to another LHE file using the `WriteToLHE` function. Using an example filename `Showeredeejj_ECM206.lhe` and the showered events in the previous example, the proper synatx is

```julia
sigma = 1.2
sstdev = 0.6
WriteToLHE(showeredEvents, "Showeredeejj_ECM206.lhe", energy, sigma, sstdev)
```

The `sigma` value is the evaluated cross section of the events and the `sstdev` is the error of it.  