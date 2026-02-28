#Permet de faire définir les courbes de Béziers en fonction de leurs points de contrôle

# POSSIBLE UPADATE
# Les noms de fonctions appelés dans intersection line doivent être plus génériques si je veux changer de paramétrisation

J = [ 0 -1. ; 1. 0]

#VALIDÉ
function bernstein(t::Float64,j::Int64,m::Int64)::Float64
# Calcule les polynomes de bernstein B^m_j(t)
# Convention B^m_j = 0 if j < 0  or j > m
    if j < 0 || j > m
        return 0.
    else
        u = t^j
        v = (1. - t)^(m-j)
        return binomial(m,j) * u * v
    end
end


#VALIDÉ
function bezier(control::Matrix{Float64},t::Float64)::Vector{Float64}
#Compute a beziers curve at t between 0 and 1 given the n control points
    n = size(control)[2] -1
    return sum( bernstein(t,k,n) * control[:,k+1] for k = 0:n )
end



# calcul le vecteur tangent à une courbe de Béziers pour un certain paramètre t
function firstderivativebezier(control::Matrix{Float64},t::Float64)::Vector{Float64}
#Compute the 'time' derivative of the beziers curve at time $t$ given the n control points
    n = size(control)[2] -1
    return sum( n*( bernstein(t,k-1,n-1) - bernstein(t,k,n-1) ) * control[:,k+1] for k = 0:n)
end

# calcul le vecteur normale à une courbe de Béziers pour un certain paramètre t
function normalbezier(control::Matrix{Float64},t::Float64)::Vector{Float64}
    n = size(control)[2] -1 # Degree of the Béziers curve   
    tangent = sum( n*( bernstein(t,k-1,n-1) - bernstein(t,k,n-1) ) * control[:,k+1] for k = 0:n)
    return J * tangent / norm(tangent)
end

#checked and deprecated
function secondderivativebezier(control::Matrix{Float64},t::Float64)::Vector{Float64}
#Deprecated
#Compute the second 'time' derivative of the beziers curve at time $t$ given the n control points
    n = size(control)[2] -1 # Degree of the Béziers curve   
    return sum( n* (n-1) * ( bernstein(t,k-2,n-2) - 2 * bernstein(t,k-1,n-2)  + bernstein(t,k,n-2) ) * control[:,k+1] for k = 0:n)
end


#calcule les quantités géométriques en un point t de la courbe de béziers c'est assez pratique
# VALIDÉ 
function geometricquantities(control::Matrix{Float64}, t::Float64)::Dict
    # utangent signie tangent Unitaire (normalisée donc)
    #Compute the  tangent vector, the  normal (as a left curve), the curvature 
    #Can be interesting to precompute the bernstein coefficients 
    n = size(control)[2] - 1  # Degree of the Béziers curve   
    point = [0. ; 0.]
    tangent = [0. ; 0.] 
    der2tangent = [0. ; 0.]
    for k in 0:n
        # Compute current point
        b = bernstein(t,k,n)
        point += b * control[:,k+1]
        # Compute tangent vector
        bp = n * ( bernstein(t,k-1,n-1) - bernstein(t,k,n-1) )
        tangent += bp * control[:,k+1]
        # Calcul second derivative 
        bpp = n * (n-1) * ( bernstein(t,k-2,n-2) - 2 * bernstein(t,k-1,n-2)  + bernstein(t,k,n-2) )
        der2tangent += bpp * control[:,k+1]
    end
    normt = norm(tangent)
    normal = J * tangent / normt #unit normal
    utangent = tangent / normt
    curvature = ( tangent[1] * der2tangent[2] - tangent[2] * der2tangent[1] ) / normt^3 # algebraic curvature at t
    return Dict("point"=> point, "utangent"=>utangent, "normal"=> normal, "curvature" => curvature, "normtangent" => normt)
end

function shapederivativeterms(control::Matrix{Float64}, t::Float64)::Dict
    n::Int64 = size(control)[2] - 1  # Degree of the Béziers curve  
    thetan = zeros(2,n+1)
    thetanprime = zeros(2,n+1)
    dic = geometricquantities(control,t)
    for j in 0:n
        thetan[:,j+1] = bernstein(t,j,n) .* dic["normal"]
        thetanprime[:,j+1] =  -dic["curvature"] .*  bernstein(t,j,n)  .* dic["utangent"]
        thetanprime[:,j+1] +=  n .* dic["normtangent"]^(-1) .* ( bernstein(t,j-1,n-1) - bernstein(t,j,n-1) ) .* dic["normal"]
    end
    return Dict("thetan" => thetan, "thetanprime" => thetanprime)
end

function basiscontrol(n::Int64, l::Int64, c::Int64)::Matrix{Float64}
    e::Matrix{Float64} = zeros(2,n)
    e[l,c] = 1.
    return e
end