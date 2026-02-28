nmale = [0., -1.] # normale au plan foca

struct Ray
    pos::Vector{Float64} #position du rayon
    angle::Float64 #direction du rayon
end

#Peut-on broadcatsre sur intersection line et lui passer un ensemble de rayons plutôt ?
# cahnger le nom line -> rayon
#On peutc changer et ne plu faire tourner les rayons
function intersectionLine(control::Matrix{Float64}, ray::Ray)::Float64
    # Les fonctions dans optim newton doivent pouvoir prendre des paramètres en compte
    n = size(control)[2] - 1 # Degree of the Béziers curve 
    theta = ray.angle
    rotation = [ cos(theta) sin(theta) ; -sin(theta) cos(theta) ] #rotation to make the line passing through x with direction xi coincides with x axis
    control0 = rotation * (control - repeat(ray.pos,1,n+1))
    function firstcomponentbeziers(t::Float64)::Float64
        return bezier(control0,t)[2]
    end
    function derfirstcomponentbeziers(t::Float64)::Float64
        return firstderivativebezier(control0,t)[2]
    end
    # approximate initial time
    t::Float64 = OptimNewton(firstcomponentbeziers, derfirstcomponentbeziers, 0.5)
    return t
end

# Il faudrait pouvoir broadcaster ausi
function billiard(control::Matrix{Float64}, ray::Ray)::Ray
    t = intersectionLine(control, ray)
    dic = geometricquantities(control, t)
    normal = get(dic,"normal", 4) 
    theta = ray.angle
    dir = [cos(theta), sin(theta)] - 2. * dot([cos(theta),sin(theta)],normal) .* normal
    thetar = atan( dir[2] , dir[1])
    rayr = Ray( bezier(control, t ), thetar)
    return rayr
end


function projectionplanfocal(ray::Ray)::Ray
# Takes a point $x$ and a direction $xi$ and returns the coordinates of impact point on the focal plane
    return Ray( [ray.pos[1] - ray.pos[2] * cot(ray.angle) , 0. ] , ray.angle)
end
#ray.pos[1] + t * cos(angle) , ray.pos[2] + t * sin(angle) 

#Differentiel de la projection sur la plan focal. Permet de définir le gradiet sur le miroir.
function diffprojectionplanfocalx(ray::Ray, dx::Vector{Float64}, dxi::Vector{Float64})::Vector{Float64}
    return dx .+  (dot([cos(ray.angle),sin(ray.angle)] , dx) ./ sin(ray.angle) ) .* nmale
end

function diffprojectionplanfocalxi(ray::Ray, dx::Vector{Float64}, dxi::Vector{Float64})::Vector{Float64}
    rayp = projectionplanfocal(ray::Ray)
    tau = norm(rayp.pos - ray.pos)
    return tau .* (dx .+  (dot([cos(ray.angle),sin(ray.angle)] , dx) ./ sin(ray.angle) ) .* nmale)  .+ dxi
end


function jacimpactpoint(d::Float64 ,  x::Vector{Float64}, xi::Vector{Float64}, dX::Vector{Float64} , dXI::Vector{Float64})::Vector{Float64}
    return dX .- (dot(nmale, dX) / dot(nmale, xi)).* xi + d .* ( dXI .- (dot(nmale, dXI) / dot(nmale, xi)).* xi )
end