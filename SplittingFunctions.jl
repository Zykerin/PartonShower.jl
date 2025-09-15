include("Constants.jl")



# The g -> gg splitting function, overestimate, integral, and integral inverse
function Pgg(z::Float64)
    return  Ca * ((1 - z * (1-z))^2) / (z * (1-z))
end
function Pgg_over(z::Float64)
    return Ca * (1/(1-z) + 1/z)
end
function Intgg(z::Float64, aSover::Float64)
    return -Ca * (aSover / (2 * pi)) * log(1/z - 1)
end
function InvInt_gg(z::Float64, aSover::Float64)
    return 1 / (1 + exp(-z/ (Ca * aSover / (2 * pi))))
end


# The g -> qqbar splitting function, overestimate, integral, and integral inverse
function Pgq(z::Float64)
    return Tr *(1- 2 * z * (1-z))
end
function Pgq_over(z::Float64)
    return Tr
end
function Intgq(z::Float64, aSover::Float64)
    return Tr * (aSover /(2 * pi)) * z
end
function InvInt_gq(z::Float64, aSover::Float64)
    return z * 2 * pi / (Tr * aSover)
end


#= The q -> gq splitting function, overestimate, integral, and integral inverse
function Pqg(z::Float64)
    return Pqq(1-z)
end
function Pqg_over(z::Float64)
    return Pqq_over(z)
end
function Int_qg(z::Float64, aSover::Float64)
    return 2 * Cf * (aSover / (2 * pi)) * log(z)
end
function InvInt_qg(z, aSover)
    return exp(2 * pi *z / (2 * Cf * aSover))
end
=#

# The q -> qg Splitting Function, overestimate, integral, and integral inverse
function Pqq(z::Float64)
    return Cf * (1 + z^2)/(1 - z)
end
function Pqq_over(z::Float64)
    return Cf * 2 / (1 - z)
end
function Intqq(z::Float64, aSover::Float64)
    return -2 * Cf * (aSover / (2 * pi)) * log(1-z)
end
function InvInt_qq(z::Float64, aSover::Float64)
    return 1- exp(-z / (2 * Cf * aSover / (2 * pi)))
end
