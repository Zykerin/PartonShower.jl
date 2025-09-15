
# The structure for an emission
struct Emission
   t::Float64
   z::Float64
   pTsq::Float64
   phi::Float64
   Generated::Bool
   continueEvolution::Bool
end

# The structure for a particle.
mutable struct Particle
    id::Int
    status::Int
    t::Float64
    z::Float64
    m::Float64
    pT::Float64
    px::Float64
    py::Float64
    pz::Float64
    E::Float64
    phi::Float64
    qT::Vector{Float64}
    color::Int32
    antiColor::Int32
    alpha::Float64
    beta::Float64
    virtuality::Float64
    continueEvolution::Bool
    aorb::String
    children::Vector{Particle}


end

struct Jet
   AllParticles::Vector{Particle}
   Progenitor::Particle
end


struct Event
   AllParticles::Vector{Particle}
   Jets::Vector{Any}


end


