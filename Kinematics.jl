include("Constants.jl")

using Roots
using LinearAlgebra

# Get the rotation matrix for two given vectors
function getRotationMatrix(v1::Vector{Float64}, v2::Vector{Float64})
    # The the cross product of two vectors
    k1 = cross(v1, v2)

    if norm(k1) < 1E-12
        return [1 0 0; 0 1 0; 0 0 1]
    else
        # Get unit vector of k
        k = k1 / norm(k1)
        # Get the k matrix
        K = [0 -k[3] k[2]; k[3] 0 -k[1]; -k[2] k[1] 0]
        
        theta = acos(clamp(dot(v1/norm(v1), v2/norm(v2)), -1, 1))

        # Get the rotation matrix
        rotMat = [1 0 0; 0 1 0; 0 0 1] + sin(theta) * K + (1 - cos(theta)) * K^2
        return rotMat
    end
end

# Function to rotate a given particle with a given rotation matrix
function rotate(p::Particle, rotMat::Matrix)
    pvec = [p.px, p.py, p.pz]
    rotVec = rotMat * pvec
    p.px = rotVec[1]
    p.py = rotVec[2]
    p.pz = rotVec[3]
end


# Function to rotate the particle to allign with the mother particle's momentum in the lab frame
function rotateMomentaLab(p::Particle, particles::Vector{Particle})
    # If the progenitor has not been radiated, then there is no reason to rotate it.
    if p.status ==1
        return
    end
    pmag = sqrt(p.px^2 + p.py^2 + p.pz^2)

    m = getRotationMatrix([0, 0, pmag], [p.px, p.py, p.pz])
    for p2 in particles
        v = [p2.px, p2.py, p2.pz]
        rotatedVec = m * v
        p2.px = rotatedVec[1]
        p2.py = rotatedVec[2]
        p2.pz = rotatedVec[3]


    end
end

# Function to get the boost factor for the outgoing jet (new) and parent jet (old)
function boostFactor(k::Float64, new::Vector{Float64}, old::Vector{Float64})
    
    qs = new[1]^2 + new[2]^2 + new[3]^2
    q = sqrt(qs)
    Q2 = new[4]^2 - new[1]^2 - new[2]^2 - new[3]^2
    kp = k * sqrt(old[1]^2 + old[2]^2 + old[3]^2)
    kps = kp^2
    betamag = (q * new[4] - kp * sqrt(kps + Q2)) / (kps + qs + Q2)

    beta = betamag * (k / kp) * Float64[old[1], old[2], old[3]]

    if betamag >= -1E-5
        return beta
    else
        return Float64[0, 0, 0]
    end

end


function boost(p::Particle, beta::Vector{Float64})
    if beta == Float64[0, 0, 0]
        return
    end

    bmag = sqrt(beta[1]^2 + beta[2]^2 + beta[3]^2)
    gamma = 1/ sqrt(1 - bmag^2)

    # Use the matrix of a lorentz boost from https://www.physicsforums.com/threads/general-matrix-representation-of-lorentz-boost.695941/
    oldp = deepcopy(p)

    p.px = - gamma * beta[1] * oldp.E + (     1 + (gamma -1) * beta[1]^2 / bmag^2) * oldp.px + ((gamma - 1) * beta[1] * beta[2] / bmag^2) * oldp.py +  ((gamma - 1) * beta[1] * beta[3] /bmag^2) * oldp.pz
    p.py = - gamma * beta[2] * oldp.E + ( (gamma -1) * beta[1] * beta[2] / bmag^2) * oldp.px + (    1 + (gamma - 1) * beta[2]^2 / bmag^2) * oldp.py + ((gamma - 1) * beta[2] * beta[3] / bmag^2) * oldp.pz
    p.pz = - gamma * beta[3] * oldp.E + ( (gamma -1) * beta[1] * beta[3] / bmag^2) * oldp.px + ((gamma - 1) * beta[3] * beta[2] / bmag^2) * oldp.py + (     1 + (gamma - 1) * beta[3]^2/ bmag^2) * oldp.pz
    p.E =    gamma * oldp.E - gamma * beta[1] * oldp.px - gamma * beta[2] * oldp.py - gamma * beta[3] * oldp.pz
end

# The function to numerically solve to find k
function kEq(k::Float64, p::Vector{Float64}, q::Vector{Float64}, s::Float64)
    sump = 0
    for i in range(1, length(p))

        sump = sump + sqrt(k^2 * p[i] + q[i])
        
    end
    sump = sump - s
    return sump
end

# Function to numerically solve for k
function solvekFactor(pj::Vector{Float64}, qj::Vector{Float64}, s::Float64)
    # Turn the function into one that can be parsed with the arguements into the root finder
    kEQ = (k -> kEq(k, pj, qj, s))
    sol = find_zero(kEQ, 0.99)
    return sol
end


# Function to perform global momemntum conservation
function globalMomCons(showeredParticles::Vector{Particle}, Jets::Vector{Jet})
    pj = Float64[] # The momenta of the parent parton
    qj = Float64[] # The momenta of the jets
    newqs = Vector{Float64}[] # Array to hold the outgoing jet's momentum
    oldps = Vector{Float64}[] # The progenitor's momentum
    rotms = Matrix{Float64}[] # The rotation matrices
    # Initialize the total energy
    sqrts = 0


    for jet in Jets
        # Append the 3-momentum of the Jet's progentior
        append!(pj, jet.Progenitor.px^2 + jet.Progenitor.py^2 + jet.Progenitor.pz^2)

        # Get the jet progenitor's 4 momentum
        oldp = [jet.Progenitor.px, jet.Progenitor.py, jet.Progenitor.pz, jet.Progenitor.E]

        # Add this jet's progenitor's energy
        sqrts += jet.Progenitor.E

        # Get the total jet's momentum after showering
        newq = Float64[0, 0, 0, 0]
        for p in jet.AllParticles
            newq += [p.px, p.py, p.pz, p.E]
        end

        # Calculate the momentum and append it to the list of jet momentum
        qj2 = newq[4]^2 - newq[1]^2 - newq[2]^2 - newq[3]^2
        if isnan(qj2)
            qj2 = 0
        end
        rotMat = getRotationMatrix([newq[1], newq[2], newq[3]], [oldp[1], oldp[2], oldp[3]])
        push!(oldps, oldp)
        push!(newqs, newq)
        append!(qj, qj2) 
        push!(rotms, rotMat)
    end

    # Get the k factor
    k = solvekFactor(pj, qj, sqrts)
    rotatedShoweredParticles = []
    # Get the electron's and postirons and append them to the list of roated showered particles since these don't need to be rotated
    for p in showeredParticles
        if abs(p.id) == 11
            push!(rotatedShoweredParticles, p)
        end
    end

    # Check if any of the jets have been radiated, if not, do not perform the boost and append the progenitor of each jet
    ifRadiated = any(length(jet.AllParticles) > 1 for jet in Jets)
    if ifRadiated == false
        for j in Jets
            push!(rotatedShoweredParticles, j.AllParticles[1])
        end
    else
        if length(Jets[1].AllParticles) == 1 || length(Jets[2].AllParticles) == 1
        end
        for (i, jet) in enumerate(Jets)
        
            # Get the boost factor
            beta = boostFactor(k, newqs[i], oldps[i])
            # Iterate through the jet's particles and rotate and boost each one
            for p in jet.AllParticles
                rotate(p, rotms[i])
                boost(p, beta)
                push!(rotatedShoweredParticles, p)

            end
        end
    end
    return rotatedShoweredParticles
end


# Function to reconstruct the sudakov basis for the entire tree with the progenitor being the starting root 
function reconSudakovBasis(prog::Particle, progPart::Particle)
    # If the progenitor particle has not emitted/on-shell, then do not reconstruct the basis
    if prog.status == 1
        return
    end
    current = prog.children[1]
    parent = prog
    stack = [prog]

    # If the stack is empty and the currently selected particle has no children, i.e. off-shell then stop the tree search
    # This loop just calculateds the relative transverse momentum and alphas
    while length(stack) != 0 || current.status == -1

        # If the current particle is off-shell i.e. there is an emisison, then append the current particle to the stack
        # then calculate the emisison of the current particle and then select the next particle as the first child in the 
        # list of children.
        if current.status == -1
            calculatePhysicals(current, prog, progPart, parent)
            push!(stack, current)
            parent = current
            current = current.children[1]
        
        # If the current particle is on-shell, i.e. no emission, then calculate its phyiscals, then select the next particle
        # as the second child of its parent and remove this parent of the stack of particles
        elseif current.status == 1
            calculatePhysicals(current, prog, progPart, parent)
            parent = pop!(stack)
            current = parent.children[2]
        end

    end
    # Algorithm misses most right particle in the tree since it ends after the current particle is set to the most right hand
    # so this calculates the basis for this particle.
    calculatePhysicals(current, prog, progPart, parent)


    current = prog
    stack = []

    # This goes through the list and calculates each beta from the final state partons up to the progenitor
    # If the stack is empty and the current off shell particle has no children whos betas are zero, we are done
    while length(stack) != 0 || (current.status == -1 && (current.children[1].beta == 0  || current.children[2].beta == 0 ))

        # Cases for when the current particle is off shell
        if current.status == -1 
            # If its first child has a zero beta then append itself to the stack and go to the first child
            if current.children[1].beta == 0 
                #print("Hit child 1 \n")
                push!(stack, current)
                current = current.children[1]
            # Check same with the second child
            elseif current.children[2].beta == 0
                #print("Hit child 2 \n")
                push!(stack, current)
                current = current.children[2]
            # If both childre have a nonzero beta, then caclulate this particles beta and go to its parent/ top in stack
            else
                calcBetai(current, prog, progPart)
                finalizeSudakov(current, prog, progPart)
                current = pop!(stack)
            end

        # If the current particle is on-shell, then claculate its beta and go to its parent particle
        elseif current.status == 1
            calcBetai(current, prog, progPart)
            finalizeSudakov(current, prog, progPart)
            current = pop!(stack)
        
        end

    end



end

# Function to calculate the beta of a given particle
function calcBetai(part::Particle, prog::Particle, progPart::Particle)
    # Case for when the particle is on-shell
    if part == prog
        
    elseif part.status == 1
        pmag = sqrt(prog.px^2 + prog.py^2 + prog.pz^2)
        nmag = sqrt(progPart.px^2 + progPart.py^2 + progPart.pz^2)

        pdotn = dot4Vec([0, 0, pmag, pmag], [0, 0, -nmag, nmag])

        betai = (part.virtuality - part.alpha^2 * part.m^2 - (-part.qT[1]^2 - part.qT[2]^2)) / (2 * part.alpha * pdotn)
        part.beta = betai
    # Case for when the particle is off-shell which is caculated iteratively from its children
    elseif part.status == -1
        betai = part.children[1].beta + part.children[2].beta
        part.beta = betai
    end
end

# Function to finally reconstruct the sudakov basis
function finalizeSudakov(part::Particle, prog::Particle, progPart::Particle)
    if part == prog
        return 
    end

    pmag = sqrt(prog.px^2 + prog.py^2 + prog.pz^2)
    nmag = sqrt(progPart.px^2 + progPart.py^2 + progPart.pz^2)
    part.px =  part.qT[1]
    part.py = part.qT[2]
    part.pz = part.alpha * pmag + part.beta * (-nmag)
    part.E = part.alpha * pmag + part.beta * (nmag)
end

# Function to actually calculate the phyiscals of the sudakov basis for a given particle
function calculatePhysicals(part::Particle, prog::Particle, progPart::Particle, parent::Particle)

    alpha = parent.alpha * part.z
    part.alpha = alpha
    
  

    if part.aorb == "c"
        kT = [-part.pT * cos(part.phi), -part.pT * sin(part.phi), 0 , 0] # pT = (px, py, pz, E)
    elseif part.aorb == "b"
        kT = [part.pT * cos(part.phi), part.pT * sin(part.phi), 0 , 0] # pT = (px, py, pz, E)
    end
    part.qT = [parent.qT[1] * part.z + kT[1], parent.qT[2] * part.z + kT[2], 0, 0] # qT = (px, py, pz, E)
    

  

end

function dot4Vec(v1::Vector{Float64}, v2::Vector{Float64})
    return v1[4] * v2[4] - v1[1] * v2[1] - v1[2] * v2[2] - v1[3] * v2[3]
end

# Function to find the color partner of a particle
function findColorPartner(parti::Particle, particles::Vector{Particle})

    partner = 0

    for pc in particles
        if parti.color == pc.antiColor && parti.antiColor == pc.color
            partner = pc
        end
    end 
    return partner

end

# Function to check that momentum is conserved for each individual particle
function checkMomCons(p::Particle)
    energy = sqrt(p.px^2 + p.py^2 + p.pz^2)

    if abs(energy - p.E) > 1E-5
        print("Momentum not conserved: momentum = " * string(energy) * ", actual = " * string(p.E) * ", difference= " * string(abs(energy-p.E)) * "\n")
    end

end

# Function to check the global momentum conservation
function checkGlobalMomCons(ev::Event)
    tot = [0, 0, 0, 0]
    for p in ev.AllParticles
        if p.status == 1
            tot += [p.px, p.py, p.pz, p.E]
        end
    end

    print("The total momentum is " * string(tot)* "\n")
end