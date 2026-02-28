# Faire des jolis plots

#POSSIBLE UPDATES
# vectoriser plot
# mettre des asserts
# mettre les plots by field !
#afiicher les points d'impact


#ne peut pas afficher qu'un seul rayon ce qui est pénible

using Plots 

function Plots.plot!(ray::Ray; choicecolor= :blue, choicelabel = false, dashreflected = false)
    rayr = billiard(control,ray)
    rayf = projectionplanfocal(rayr)
    timpact1 = norm(rayr.pos-ray.pos)
    timpact2 = norm(rayr.pos - rayf.pos) 
    theta = ray.angle
    thetar = rayr.angle
    x = ray.pos[1] -0.1 * timpact1 * cos(theta)
    xr = rayr.pos[1]
    y = ray.pos[2] -0.1 * timpact1 * sin(theta)
    yr = rayr.pos[2] 
    xf = rayf.pos[1]
    yf = rayf.pos[2]
    if dashreflected == true
        # Plot du rayon incident
        plot!([x , xr ], [y , yr]; color = choicecolor, label = choicelabel, linestyle = :dash)
        # Plot du rayon réflechi
        plot!([xr , xf ], [yr, yf]; color = choicecolor, label = choicelabel)
    else
        plot!([x , xr, xr , xf ], [y, yr, yr, yf]; color = choicecolor, label = choicelabel)
    end
end

function circleshape(xcenter::Float64,ycenter::Float64,r::Float64)::Tuple{Vector{Float64},Vector{Float64}}
    theta = [ 2. * pi * i / 500. for i in 0:500]
    X = xcenter .+ r .* cos.(theta)
    Y = ycenter .+ r .* sin.(theta)
    return (X,Y)
end

function segment(angle::Float64, centre::Float64, longueur::Float64)::Tuple{Vector{Float64},Vector{Float64}}
    Y = [ centre - longueur + 2 * i * longueur / 500. for i in 0:500]
    X =  0. .* Y .+ angle
    return (X,Y)
end
