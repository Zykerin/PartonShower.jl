
# The QCD constants
global const Nc = 3 # The number of colors
global const Cf = (Nc^2 -1)/ (2 * Nc) # The quark color factor associated with gluon emissions from a quark
global const Ca= Nc # The color factor associated with gluon emissions from another gluon
global const Tr = 1/2 # The color facor for a gluon to qqbar emission



function initializeAlphaS(asmz0, mz0, mb0=4.2, mc0 =1.25, mt0= 174.2, order0=2)
    global asmz = asmz0
    global mz = mz0
    global mb = mb0
    global mc = mc0
    global mt = mt0
    global order = order0
    global lambdas = computeLambdas()


end

function beta0(nf::Int64)
    return (11 / 6 * Ca) - (2 / 3 * Tr * nf)
end

function beta1(nf::Int64)
    return (17 / 6 * Ca^2) - ((5 / 3 * Ca + Cf) * Tr * nf)
end


function lambdaFromAs(Q::Float64, asQ::Float64, nf::Int64)
    b0 = beta0(nf) / (2 * pi)
    b1 = beta1(nf) / (2 * pi)^2
    t = 1 / (b0 * asQ)
    if order == 1
        return Q * exp(-t/2)
    elseif order == 2
        return Q * exp(-t/2) * (b0 * asQ)^(-b1/ (2 * b0^2))
    else 
        throw("Order must be 1 or 2.")
    end
end

function asFromLambda(Q::Float64, lambda::Float64, nf::Int64)
    b0 = beta0(nf) / (2 * pi)
    b1 = beta1(nf) / (2 * pi)^2
    t = log(Q^2 / lambda^2)
    if order == 1
        return 1 / (b0 * t)
    elseif order == 2
        return 1 / (b0 * t) * (1 - b1 * log(t)/ (b0^2 * t))
    else
        throw("Order must be 1 or 2.")
    end
end


function computeLambdas()
    lambdas2 = Dict()
    # Step 1: Lamba5 from mz
    lambda5 = lambdaFromAs(mz, asmz, 5)
    lambdas2[5] = lambda5
    # Step 2: Match down to nf=4 at mb
    asmb5 = asFromLambda(mb, lambda5, 5)
    lambda4 = lambdaFromAs(mb, asmb5, 4)
    lambdas2[4] = lambda4
    # Step 3: Match down to nf=3 at mc
    asmc4 = asFromLambda(mc, lambda4, 4)
    lambda3 = lambdaFromAs(mc, asmc4, 3)
    lambdas2[3] = lambda3
    # Step 4: Match up to nf = 6 at mt (if needed)
    asmt5 = asFromLambda(mt, lambda5, 5)
    lambda6 = lambdaFromAs(mt, asmt5, 6)
    lambdas2[6] = lambda6
    return lambdas2
end



function alphaQ(Q::Float64)
    if Q < mc
        nf = 3
    elseif Q < mb
        nf = 4
    elseif Q < mt
        nf = 5
    else
        nf = 6
    end
    lambda = lambdas[nf]
    return asFromLambda(Q, lambda, nf)
    
end




function getAlphaS(t::Float64, z::Float64, Qcut::Float64)
    scale = z * (1-z) * sqrt(t)
    if scale < Qcut
        scale = Qcut
    end
    return alphaQ(scale)
end