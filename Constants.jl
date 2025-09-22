include("alphaS_HW.jl")

global const Qc = 0.935

# Initialize the alphaS variable
initializeAlphaS(0.1074, 91.1876)

aSover::Float64 = alphaQ(Qc) # The overestimate for the coupling constants

global const pTmin::Float64 = 0.900
global const pT2min::Float64 = pTmin^2


a = 0.3
b = 2.3
delta = 2.8
c = 0.3

# These are the functions used for the masses of the quarks and gluons
function Qg(mq::Float64)
   return max((delta - a * mq)/b, c)

end

function mufunc(mass::Float64, ps::Float64)
   return max(mass, Qg(ps))
end
