# Fonctions d'optimisation

######### NEWTON' METHOD ########
# look for the zero of func using Newton method

# funcder is the derivative of func
# init is the initial point where we compute the tangent
# numberier is a number of iterations not to exceed
# precision is a treshold number carrying on func(x) below which we stop the computation

struct NoSolutionException <: Exception end

# POSSIBLE UPDATES
# Mettre le gradient projeté (modification mineur du gradient)
#Optim Newton ne marche que pour des fonctions de la variable réel pour le moment
function OptimNewton(func::Function, funcder::Function, x0::Float64 , numberiter::Int64 = 30, precision::Float64 = 10^(-15) , verbose::Bool = true)::Float64
    j::Int64 = 0
    xold = x0 + 1.
    while abs(x0 - xold) > precision  && j < numberiter
        # can put another criterion on distance between points between iteration (slow convergence)
        j += 1
        if verbose == true
            #println("Iteration number : ", j)
            #println("Value of f(x) at the current iteration : ", func(x0))
        end
        xold = x0
        x0 = x0 - func(x0) / funcder(x0)
    end
    if (j == numberiter) && (verbose=true)
        println("Max iterations reached")
        println("Newton may not have solutions")
        throw(NoSolutionException)
    elseif verbose == true
        #println("Newton has converged at iteration " , j)
        #println("The value of f is " , func(x0))
    end
    return x0
end

######### GRADIENT METHOD ########
# return a point x where func is minimized (or inflexion point)

#Mettre Array pour l'indiquer à OptimGradient Array est un cas général
function OptimGradient(func::Function, funcder::Function, x0::Array, step::Float64, distance1::Function, distance2::Function; proj::Function = x -> x, history::Bool = false, numberiter::Int64 = 100, iterstepsearch::Int64 = 15 ,  verbose::Bool = false)::Array
    # distance1 is the distance to measure how far the gradient is from 0 typically a norm on vectors of R^d
    # distance2 is the distance to measure the distance between to consecutive steps typically the absolute value
    # proj est utile pour minimiser sur un convexe fermé
    f = func(x0)
    if verbose
        println("The initial value of the cost is : ", f)
    end
    h = step # h is the step to update in internal descent
    
    df =  funcder(x0) # derivative of func    
    x = proj(x0 .- h .* df) # x is computed iteratively
    
    k = 0 # stopping index for algorithm too long to converge
    veps1 = 10^(-16) #critère gradient
    veps2 = 10^(-10)
    c = 0.5 # Armijo
    if history
        # iter et cost sont des variables globales
        append!(iter, 0)
        append!(cost,f)
    end
    while distance1(df) >= veps1 && k < numberiter 
        if verbose
            println("Iteration number : ", k)
            println("Value of x : ", x)
            println("Value of nabla f(x) at the current iteration : ", funcder(x0))
        end
        k += 1
        j = 0
        
        while func(x) >= f - c * h * norm(df)^2 && j < iterstepsearch # Update of x Armijo
            j += 1
            h = h ./ 2 
            x = proj(x0 .- h .* df)  # descent : should work in theory
        end
        # Pour sauvergarder l'histoire de convergence
        if history
            append!(iter, k)
            append!(cost,func(x))
        end
        #If max iterations reached
        if j == iterstepsearch
                println("Warning : Algorithm may not have converged")
                println(" Number of iterations for descent research exceeded !")
                println("The value of the gradient norm (distance 1) before stopping is : ", distance1(funcder(x)))
                println("The value of the cost before stopping is : ", func(x))
            return x
        end
        if distance2(x0 .- x) < veps2 
                println("Warning : Algorithm may not have converged")
                println("Distance between two consecutive updates is small")
                println("The final value of the gradient norm (distance 1) is : ", distance1(funcder(x)))
                println("The final value of the cost is : ", func(x))
            return x
        end
        x0 = x 
        f = func(x0)
        df = funcder(x0)
        h = 1.7 * step # we let ourselves a chance to increase the steps 
        x = proj(x0 .- h .* df) # x is computed iteratively

    end
    if k == numberiter
        println("Warning : Algorithm may not have converged and you can increase the number of total steps")
        println("Number of total iterations exceeded")
        println("The final value of the gradient norm (distance 1) is : ", distance1(funcder(x)))
        println("The final value of the cost is : ", func(x))
        return x
    elseif distance1(df) < veps1
        if verbose
            println("User convergence treshold attained ! ")
            println("The value of derfunc for the returned x is :", distance1(funcder(x)))
            println("The value of func for the returned x is :", func(x) )
        end
        return x
    end
end


# on reshape df avant d'inverser puis on déreshape alpha après l'inversion sinon on garde tout pareil "2 * 9"

function LMdescent(func::Function, funcder::Function,  xold::Array,  distance::Function,  proj::Function = x -> x, step0::Float64 =0.1,  history::Bool = false, numberiter::Int64 = 100, iterstepsearch::Int64 = 15 ,  verbose::Bool = false)::Array
    # distance is the distance to measure the distance between to consecutive steps typically the absolute value
    l::Int64 = size(xold)[2]
    k::Int64 = 0 # stopping index for algorithm too long to converge
    veps::Float64 = 1e-10
    if history
        # iter et cost sont des variables globales
        f = func(xold)
        append!(iter, 0)
        append!(cost,f)
    end
    xold = vcat(xold[1,:] , xold[2,:])
    x = xold .+ 1.
    lambda = 1. / step0
    while distance(xold .- x)  >= veps && k < numberiter
        println("Iteration number : ", k)
        if k >= 1
            xold = x
        end
        k += 1
        j = 0
        f = func(Matrix(reshape(xold,l,2)'))
        if verbose
            println("Iteration number : ", k)
            println("Value of x : ", x)
            println("Value of nabla f(x) at the current iteration : ", funcder( Matrix(reshape(x,l,2)')) )
        end
        (df,dr2) = funcder(Matrix(reshape(xold,l,2)'))
        df = proj(df)
        gf = vcat(df[1,:] , df[2,:])
        alpha = ( dr2 +  UniformScaling(lambda) ) \ gf
        x = xold .- alpha
        while func( Matrix(reshape(x,l,2)') ) >= f  && j < iterstepsearch 
            j += 1
            lambda = 2. * lambda
            alpha = ( dr2 +  UniformScaling(lambda)) \ gf
            x = xold .- alpha  
        end
        lambda = lambda / 5.
        println("Le lambda est :",lambda)
        println("La plus petite valeur propre est :", eigmin(dr2))
        # Pour sauvergarder l'histoire de convergence
        if history
            append!(iter, k)
            append!(cost,func(Matrix(reshape(x,l,2)')))
        end
        #If max iterations reached
        if j == iterstepsearch
                println("Warning : Algorithm may not have converged")
                println(" Number of iterations for descent research exceeded !")
                println("The value of the cost before stopping is : ", func(Matrix(reshape(x,l,2)')) )
            return Matrix(reshape(x,l,2)')
        end
    end
    if k == numberiter
        println("Warning : Algorithm may not have converged and you can increase the number of total steps")
        println("Number of total iterations exceeded")
        println("The final value of the cost is : ", func( Matrix(reshape(x,l,2)') ) )
        return Matrix( reshape(x,l,2)' )
    end
    if verbose
        println("User convergence treshold attained ! ")
        println("The value of func for the returned x is :", func( Matrix(reshape(x,l,2)') ) )
    end
    println("Tolerance on the distance between two consecutive iterates reached :", distance(xold .- x) )
    return Matrix( reshape(x,l,2)' )
end



