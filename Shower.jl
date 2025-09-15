include("ShowerHelpers.jl")
include("SplittingFunctions.jl")
include("Kinematics.jl")
include("Constants.jl")

using Random 


Random.seed!(1234) # Set the seed for easier testing

function generateEmissions(Q::Float64, Qcut::Float64, aSover::Float64, branchType::Int8, masses::Vector{Float64})

    R1 = rand()
    R2 = rand()
    R3 = rand()
    R4 = rand()
    

    t = 0
    z = 0
    Generated = true
    continueEvolution = true
    pTsq = 0
    phi = 0
    # Check for the branch type and get the appropiate z and t emissions
    if branchType == 1
        tmin = 3.99 * (masses[1]^2 + Qcut^2)
        t, continueEvolution = tEmission(Q, R1, aSover, tmin, Intgg, masses, branchType)
        z = zEmission(t, R2, aSover, Intgg, InvInt_gg, masses, branchType)
    elseif branchType == 2
        t0 = masses[1]^2 + Qcut^2
        t1 = masses[2]^2 + Qcut^2
        tmin = 0.99 * (2 * sqrt(t1 * t0) + t1 + t0)
        t, continueEvolution = tEmission(Q, R1, aSover, tmin, Intqq, masses, branchType)
        z = zEmission(t, R2, aSover, Intqq, InvInt_qq, masses, branchType)
    elseif branchType == 3
        tmin = 3.99 * (masses[1]^2 + Qcut^2)
        t, continueEvolution = tEmission(Q, R1, aSover, tmin, Intgq, masses, branchType)
        z = zEmission(t, R2, aSover, Intgq, InvInt_gq, masses, branchType)
    end
    
    # If the evoltion has been terminated, return the emission struc
    if continueEvolution == false
        z = 1
        return Emission(t, z, pTsq, phi, Generated, continueEvolution)
    end

    # Get the transverse momentum squared
    pTsq = transversemmSquared(t, z, branchType, masses)

    # Apply the minimum transverse momentum cut
    if pTsq < pT2min
        Generated = false
    end

    # Apply the overestimate of the splitting function cut according to the branch
    if branchType == 1
       if R3 > Pgg(z) / Pgg_over(z)
            Generated = false
       end
    elseif branchType == 2
        if R3 > Pqq(z) / Pqq_over(z)
            Generated = false
       end
    elseif branchType == 3
        if R3 > Pgq(z)/ Pgq_over(z)
            Generated = false
       end
    end

    # Apply the alphaS oversitmate cut
    if R4 > getAlphaS(t, z, Qcut) / aSover
        Generated = false
    end

    phi = (2 * rand() - 1) * pi

    if Generated == false
        #print("Not Generated")
        z = 1
    end

    return Emission(t, z, pTsq, phi, Generated, continueEvolution)

end

function genEmissions2(Q::Float64, Qcut::Float64, aSover::Float64, branchType::Int8, masses::Vector{Float64})
    R1 = rand()
    R2 = rand()
    R3 = rand()
    R4 = rand()
    t = 0
    generated = true
    continueEvolve = true
    z = 0
    pTsq = 0
    phi = 0
    # Check for the branch type and get the appropiate z and t emissions
    if branchType == 1
        t, continueEvolve = gettEmisisonDirect(Q, R1, Qcut, aSover, Intgg, masses, branchType)
   
    elseif branchType == 2
        t, continueEvolve = gettEmisisonDirect(Q, R1, Qcut, aSover, Intqq, masses, branchType)
     
    elseif branchType == 3
        t, continueEvolve = gettEmisisonDirect(Q, R1, Qcut, aSover, Intgq, masses, branchType)

    end
   
    if continueEvolve == false
        return Emission(t, z, pTsq, phi, generated, false)
    end
    zup, zlow = zBounds(masses, t, branchType)


    # Ensure that the limits are nonnegative
    if zup < 0 || zlow < 0
        generated = false
    end
    if branchType == 1
        z = zEmission(Q^2, R2, aSover, Intgg, InvInt_gg, masses, branchType)
    elseif branchType == 2
        z = zEmission(Q^2, R2, aSover, Intqq, InvInt_qq, masses, branchType)
    elseif branchType == 3
        z = zEmission(Q^2, R2, aSover, Intgq, InvInt_gq, masses, branchType)
    end
    # Ensure that the upper bounds is greater than the lower bounds
    if zlow > zup
        generated = false
    end

    # Check the ensure that z is within the limits
    if z > zup || z < zlow
        generated = false
    end

    # Get the transverse momentum squared
    pTsq = transversemmSquared(t, z, branchType, masses)

    # Check that the momentum is greater than the momentum cutt off
    if pTsq < pT2min
        generated = false
    end

    # Apply the overestimate of the splitting function cut according to the branch
    if branchType == 1
       if R3 > Pgg(z)/ Pgg_over(z)
            generated = false
       end
    elseif branchType == 2
        if R3 > Pqq(z) / Pqq_over(z)
            generated = false
       end
    elseif branchType == 3
        if R3 > Pgq(z)/ Pgq_over(z)
            generated = false
       end
    end

    # Apply the alphaS oversitmate cut
    if R4 > getAlphaS(t, z, Qcut) / aSover
        generated = false
    end

    phi = (2 * rand() - 1) * pi


    if generated == false
        z = 1
    end

    return Emission(t, z, pTsq, phi, generated, continueEvolve)

end

# Functionn to evolve a specific particle given its children
function evolveParticle(pa::Particle, pb::Particle, pc::Particle, Qcut::Float64, aSover::Float64)
    tcuttoff = 4 # The curtoff for the evolution variable

    emission = Emission(pa.t, 1, 0, 0, true, true) # Base values for the emisisons
    Q::Float64 = sqrt(pa.t) * pa.z # Get the initial scale for this evolution
    masses = Float64[0, 0] # Get the masses which are set to 0 for now
    branchType::Int8 = 0
    # This loops until either an emission is accepted and generated or until the evolution is terminated
    while true
        # Condition if the evolution scale is below the cuttoff to break
        if Q < sqrt(tcuttoff * Qcut^2)
            pa.continueEvolution = false
            pa.status = 1 # Set the parent particle to final state if there is not emission
            pa.virtuality = 0 # The virtuality is zero if the particle has no emission
            return
        end
        # Branch for quark emiting. Currently only q -> qg
        if abs(pa.id) < 6 && abs(pa.id) > 0
            branchType = 2
            #pb.m = mufunc(pa.m, Qg(pa.m))
            #pc.m = Qg(pa.m)
            #emission = generateEmissions(Q, Qcut, aSover, branchType, masses)
            emission = genEmissions2(Q, Qcut, aSover, branchType, masses)
            pb.id = pa.id
            pc.id = 21
        # Condition for gluon emittion which needs competition
        elseif abs(pa.id) == 21
            branchType = 1
            #pb.m = Qg(pa.m)
            #pc.m = Qg(pa.m)

            #emission = generateEmissions(Q, Qcut, aSover, branchType, masses)
            emission = genEmissions2(Q, Qcut, aSover, branchType, masses)
            pb.id = 21
            pc.id = 21

            if emission.continueEvolution == false
                emission = Emission(0, 1, 0, 0,  false, false)
            end
            # Test for the possibility of a g -> qqbar emission
            branchType  = 3
            for flavor in range(1, 5)
                #emissionTemp = generateEmissions(Q, Qcut, aSover, branchType, masses)
                emissionTemp = genEmissions2(Q, Qcut, aSover, branchType, masses)

                # Accept the emission if the generated t value is greater 
                if emissionTemp.continueEvolution == true && emissionTemp.t > emission.t
                    emission = emissionTemp
                    #pb.m = mufunc(pa.m, Qg(pa.m))
                    #pc.m = mufunc(pa.m, Qg(pa.m))
                    pb.id = flavor
                    pc.id = -flavor
                end
            end

        end
        if emission.t < (tcuttoff * Qcut^2)
            pa.continueEvolution = false
            pa.status = 1 # Set the parent particle to final state if there is not emission
            pa.virtuality = 0 # The virtuality is zero if the particle has no emission
            return
        end
        # Rescale the evolution variable
        Q = sqrt(emission.t)
        # Condition to break out of the loop: 
        # either an emission has been generated or the evolution is terminated
        if emission.Generated == true || emission.continueEvolution == false
            break
        end

    end


    if emission.continueEvolution == false
        pa.continueEvolution = false
        pa.status = 1 # Set the parent particle to final state if there is not emission
        pa.virtuality = 0 # The virtuality is zero if the particle has no emission
        return
    end
    pb.t = emission.t
    pc.t = emission.t
    pb.z = emission.z
    pc.z = 1 - emission.z 

    pb.m = 0
    pc.m = 0

    # Now apply the color structure 
    newcolor = currentcolor + 1
    global currentcolor = newcolor
    # Case for when the emitting particle is a gluon
    if abs(pa.id) == 21
        # Case for g -> gg
        if abs(pb.id) == 21

            pb.antiColor = pa.antiColor
            pb.color = newcolor

            pc.color = pa.color
            pc.antiColor = newcolor
        # Case for g -> qqbar
        else 
            pb.color = pa.color
            pc.antiColor = pa.antiColor
        end
    # Case for  quark emitting, so currently only q -> qg
    else
        # Need to check whether emitting quark is an antiquark or not
        # For antiquark, qbar -> qbarg
        if pa.id < 0
            pb.antiColor = newcolor
            pc.antiColor = pa.antiColor
            pc.color = newcolor
        # Case for q -> qg splitting
        else
            pb.color = newcolor
            pc.color = pa.color
            pc.antiColor = newcolor
        end

    end

    pb.phi = emission.phi 
    pc.phi = emission.phi

    pT = sqrt(emission.pTsq)
    vmsq = getVirtMsq(emission.t, emission.z) # Get the virtuality of the emitting particle
    pmag = sqrt(pa.px^2 + pa.py^2 + pa.pz^2)

    pa.virtuality = vmsq
    #=
    pb.px = pT * cos(emission.phi)
    pb.py = pT * sin(emission.phi)
    pb.pz = pb.z * pmag
    pb.E = sqrt(pb.px^2 + pb.py^2 + pb.pz^2)

    pc.px = - pT * cos(emission.phi)
    pc.py = - pT * sin(emission.phi)
    pc.pz = pc.z * pmag
    pc.E = sqrt(pc.px^2 + pc.py^2 + pc.pz^2)
    =#

    pa.status = -1 # Set the emitting particle Particle's status to -1 so as to indicate that it is an intermediate particle

    pb.aorb = "b"
    pc.aorb = "c"

    pb.pT = pT
    pc.pT = pT 

    pb.status = 1
    pc.status = 1
    pb.continueEvolution = true
    pc.continueEvolution = true 

    return
end


function showerParticle(jet::Jet, particle::Particle, Qmin::Float64, aSover::Float64)
    # Append the progenitor to the jet's list of particles
    push!(jet.AllParticles, particle)

    i = 1
    while i < 26
        # Check if the current index is out of bounds for the list of particles in the jet
        # If so, we are done with showering this jet
        if i > length(jet.AllParticles)
            return
        end
        pa = jet.AllParticles[i]
        # Get the child particles' templates
        pb = Particle(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, [], 0, 0, 0, 0, 0, true, "", [])
        pc = Particle(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, [], 0, 0, 0, 0, 0, true, "", [])
        evolveParticle(pa, pb, pc, Qmin, aSover)
        # If the evolution is terminated, end this branch
        if pa.continueEvolution == false
            i += 1
            continue
        end
        # If not then set the current index to the b particle and append the c particle to the list
        jet.AllParticles[i] = pb
        push!(jet.AllParticles, pc)
        # Append the child particles, b and c, to the list of children in the parent particle
        append!(pa.children, [pb, pc])

    end

end


function showerEvent(event::Event, Qmin::Float64, aSover::Float64)

    newerEvent= Event([], [])
    plist = Particle[]
    jets = Jet[]
    AllParticles = Particle[]
    oldEvent = deepcopy(event) # Make a copy of the event in case a full event veto needs to be done

    # Go through the list of particles in the event and get the ones that will be showered and append them to a list
    # This is needed since both particles are needed for the initial evolution scale
    for p in event.Jets
        # Check if the particle is an electron or positron to skip it
        if abs(p.id) == 11
            push!(AllParticles, p)
        elseif abs(p.id) < 6 && abs(p.id) > 0 && p.status == 1 # Only shower final state particles
            push!(plist, p)
        end 
    end

    # Set the initial evolution scale for the two progenitors
    plist[1].t, plist[2].t = EvolutionScale(plist[1], plist[2])

    # Get the current larges color value for this event
    global currentcolor = maximum([p.color for p in plist])

    for (i, p) in enumerate(plist)
        jet = Jet([], p)

        showerParticle(jet, p, Qmin, aSover)
        partColor = findColorPartner(p, plist)
        
        reconSudakovBasis(p, partColor)
        
        rotateMomentaLab(p, jet.AllParticles)
        
        append!(AllParticles, jet.AllParticles)
        push!(jets, jet)
    end

    # Test to see if the momentum reconstruction has worked. If an error is thrown, ususually from the finding k part, then the reconstruction has failed 
    # and the event needs to be vetoed. 
    try 
        append!(newerEvent.AllParticles, globalMomCons(AllParticles, jets))
    catch
        newerEvent = Event(oldEvent.Jets, oldEvent.Jets)
    end


    # Check for if any of the reconstructed values are NaN 
    checkNaN = []
    for p in newerEvent.AllParticles
        if isnan(p.px)
            print("px= " * string(p.px) * ", py= " * string(p.py) * ", pz= " * string(p.pz)*", E= " * string(p.E) * "\n")
            append!(checkNaN, true)
        else
            append!(checkNaN, false)
        end
    end

    # If they are, then reshower this event
    if any(checkNaN)
        print("Reshowering \n")
        newerEvent = showerEvent(oldEvent, Qmin, aSover)
    end

    return newerEvent
end
