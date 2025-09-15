include("Structures.jl")
using StatsBase

global const alpha::Float64 = 1/132.507
global const ECM::Float64 = 206
global const s::Float64 = ECM^2
# Define GeV to pb conversion
global const pb::Float64 = 3.894E8

# Define the Fermi Constant
global const Gf::Float64 = 1.16639E-5
# Define the Z Boson mass
global const MZ::Float64 = 91.188
# Define the Z Boson Width
global const gammaZ::Float64 = 2.4414
# The electric charge for partciles
function Qf(which::Int64)

    if which == 1
        return 2/3
    elseif which == 2
        return -1/3
    elseif which == 3
        return 0
    elseif which == 4
        return -1
    end

end

# Define the vector couplings
function Vf(which::Int64)

    if which == 1
        return 1/2 - (4/3) * 0.222246
    elseif which == 2
        return -1/2 + (2/3) * 0.222246
    elseif which == 3
        return 1/2
    elseif which == 4
        return -1/2 + 2 * 0.222246
    end

end

function Af(which::Int64)

    if which == 1
        return 1/2 
    elseif which == 2
        return -1/2
    elseif which == 3
        return 1/2
    elseif which == 4
        return -1/2
    end

end

global const kappa::Float64 = (sqrt(2) * Gf * (MZ^2))/(4 * pi * alpha)
global const chi1::Float64 = (kappa * s * (s - MZ^2)) / ((s - MZ^2)^2 + gammaZ^2 * MZ^2)
global const chi2::Float64 = (kappa^2 * s^2) / ((s - MZ^2)^2 + gammaZ^2 * MZ^2)


function A0(which::Int64)
    return  Qf(which)^2 - 2 * (Qf(which)) * Vf(which) * Vf(4) * chi1 + (Af(4)^2 + Vf(4)^2) * (Af(which)^2 + Vf(which)^2) * chi2
end

function A1(which::Int64)
    return -4 * Qf(which) * Af(4) * Af(which) * chi1 + 8 * Af(4) * Vf(4) * Af(which) * Vf(4) * chi2
end

# Function for the different charge of the quark
function eq(q::Int64)
    # For up, charm, and top quarks
    if q == 1
        return 2/3
    # For down, strange, and bottom quarks
    elseif q == 2
        return -1/3
    end
end


# The differential cross section function 
function dsigma(x::Float64, q::Int64)
    return 3 * 2 * pi * (alpha^2)/ (4 * s) * (A0(q) * (1 + x^2) + A1(q) * x) 
end

function f2(x::Float64, q::Int64)
    return 3 * 2 * pi * (alpha^2)/ (4 * s) * ( (1 + x^2)) * eq(q)^2
end

# The function to perform Monte Carlo integration
function mcInt(func, x1, x2, N)

    sumw = 0
    sumwsq = 0
    maxW = -1E99
    xmaxW = -1E99

    # Perform the required number of integrations
    for i in range(1, N)
        xi = (x2 - x1) * rand() + x1
        # Set the total of the sum to 0
        functot = 0
        # Sum over the quark favors 
        for i in range(1, 4)
            if i == 1
                q =1
            elseif i == 2
                q = 2
            elseif i == 3
                q =2
            elseif i == 4
                q = 1
            else
                q =2
            end
            functot += func(xi, q)
        end
        functot
        yi = (x2 - x1) * functot
        sumw += yi
        sumwsq += yi^2

        if yi > maxW
            maxW = yi 
            xmaxW = xi
        end
    end

    I = sumw / N
    error = sqrt((1/N)* ((1/N) * sumwsq - I^2))

    return I, error, maxW, xmaxW
end



#=
global const N = 1E8

# Get the integrated cross section and compare it to the analytical for accuracy
Integ, Err, maxW, xmaxW = mcInt(dsigma, -1, 1, N)
global actual = 0
for i in range(1, 4)
    if i == 1
        q =1
    elseif i == 2
        q = 2
    elseif i == 3
        q =2
    elseif i == 4
        q = 1
    else
        q =2
    end
    global actual += 4 * pi * ((alpha^2) / (s )) * A0(q)
    #global actual += 4 * pi * ((alpha^2 * kappa^2) / ( gammaZ^2 )) * (Af(q)^2 + Vf(q)^2) * A0(q)
end

print("The result from Monte Carlo integration is " * string(Integ*pb) * " +- " * string(Err* pb) * ", while the analytical is " * string(actual*pb) * ".\n" )
print("The difference is " * string(abs(Integ*pb - actual*pb)) * "pb. \n")

=#


function eventGen(func, x1, x2, nGen, nInt)
    IntAllFlav, Err, maxW, xmaxW = mcInt(dsigma, x1, x2, nInt)
    count = 0

    eventList = []
    # Keep attempting to generate events until the number of required events is met
    while count < nGen
        xi = (x2 - x1) * rand() + x1
        # Set the total of the sum to 0
        functot = 0
        # Sum over the quark favors 
        for i in range(1, 4)
            if i == 1
                q =1
            elseif i == 2
                q = 2
            elseif i == 3
                q =2
            elseif i == 4
                q = 1
            else
                q =2
            end
            functot += func(xi, q)
        end
        yi = (x2 - x1) * functot
        # Condition to accept the event
        if yi/ maxW > rand()
            append!(eventList, xi)
            count +=1
        end

    end

    return eventList, IntAllFlav, Err
end



# Function to reconstruct the momenta of the event 
function reconMomenta(events)
    reconEvents = []
    energy = sqrt(s) /2

    # Go through each event and reconstruct it in the form of a particle class
    for costh in events
        # Just need a random phi since it doesn't particlularly matter, just as long as its consistently used
        phi = rand() * 2 * pi
        sinth = sqrt(1-costh^2)
        # Make the electron and positron
        elec = Particle(11, -1, 0, 0, 0, 0, 0, 0, -energy, energy, 0, [], 0, 0, 0, 0, 0, false, "", [])
        pos = Particle(-11, -1, 0, 0, 0, 0, 0, 0, energy, energy, 0, [], 0, 0, 0, 0, 0, false, "", [])

        weights = Float64[]
        totsigma = 0
        # Go through the quark flavors and get the total sum of the differential cross section
         for i in range(1, 4)
            if i == 1
                q =1
            elseif i == 2
                q = 2
            elseif i == 3
                q =2
            elseif i == 4
                q = 1
            else
                q =2
            end
            # Get the weight for this specific quark
            append!(weights, dsigma(costh, q))
            totsigma += dsigma(costh, q)
        end
        items = [1, 2, 3, 4]
        # Choose the quark flavor
        flavor = sample(items, Weights(weights))

        # Create the q and q bar particles
        q = Particle(flavor, 1, 0, 1, 0, 0, -energy * sinth * cos(phi), -energy * sinth * sin(phi), -energy * costh, energy, phi, [0, 0, 0, 0], 501, 0, 1, 0, 0, true, "", [])
        qbar = Particle(-flavor, 1, 0, 1, 0, 0, energy * sinth * cos(phi), energy * sinth * sin(phi), energy * costh, energy, phi, [0, 0, 0, 0], 0, 501, 1, 0, 0, true, "", [])

        push!(reconEvents, Event([], [pos, elec, q, qbar]))
    end

    return reconEvents
end

# Function to actually generate the events
function generateEvents(nEvents, nInt)
    # Generate the events
    events, Integ, error = eventGen(dsigma, -1, 1, nEvents, nInt)
    # Reconstruct the momenta
    finalEvents = reconMomenta(events)

    return finalEvents, Integ*pb, error*pb
end