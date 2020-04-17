using LinearAlgebra
using SparseArrays
using MathProgBase
using Clp

"""
    ratiotest(x, A, b, d)
Perform the ratio test: given the direction `d` choose
`α = argmin b -  A(x+alpha d) ≧ 0`.
"""

function ratiotest(x, A, b, d)
    min_ratio = 0.
    Atimesd = A*d
    count = 0
    for (index,aTd) in  enumerate(Atimesd)
        if aTd > 1e-8
           count += 1
           ratio  = (b[index] - dot(A[index,:],x))/aTd
           if ratio < min_ratio || count == 1
               min_ratio = ratio
           end
       end
    end
    return min_ratio
end


"""
    barydirection(x, A, b, d)
Finds Barycentric direction

"""
function barydirection(x, A, b, Anormed,cnormed)
    # Find the Opposite Barycenteric direction
    # Find indexes J such that b-Ax < ϵ
    J = findall(abs.(b-A*x) .< eps())
    cardJp1 = length(J) + 1.
    d = cnormed
    for index in J
        d += Anormed[index,:]
    end
    d /= -cardJp1
    return d
end

"""
    finddirection(x, A, b, Anormed,cnormed)
Finds Barycenteric direction

"""
function finddirection(x, A, b, Anormed,cnormed, num_var, sortedJ)
    # Find the Opposite to Circumcenter direction
    # Finds indexes J such that b-Ax < ϵ
    J = findall(abs.(b-A*x) .< 1e-8)
    lenJ = length(J)
    # If x is interior, takes direction -c
    if isempty(J)
        return -cnormed
    end
    index = 1
    while lenJ >  (num_var-1)
        filter!(x->x != sortedJ[index],J)
        lenJ = length(J)
        index +=1
    end
    X = Matrix([cnormed (Anormed[J,:])'])
    xcirc = FindCircumcentermSet(X)
    return -xcirc
end

"""
    refinesolution(x, A, b, c, num_var)
Refine solution when near the LP solution

"""
function refinesolution(x, A, b, c, num_var, atol)
    index_active = findall(b - A*x .<  1e-8)
    num_active = length(index_active)
    iter = 0
    if iszero(num_active)
        alpha = ratiotest(x,A,b,-c)
        x = x - alpha*c
        index_active = findall(b - A*x.<= 1e-8)
        num_active = length(index_active)
    end
    while num_active < num_var && iter <= num_var
       iter += 1
       aFact = cholesky(A[index_active,:]*A[index_active,:]')
       lambda = aFact\(A[index_active,:]*c)
       d = -c +  A[index_active,:]'*lambda
       # if norm(d) ≈ 0
       #     break
       # end
       alpha = ratiotest(x,A,b,d)
       xnew = x + alpha*d
       if norm(xnew - x) < atol
           return xnew
       end
       x = xnew
       index_active = findall(b - A*x .<  1e-8)
       num_active = length(index_active)
    end
    return x
end


"""
mpstomatrix(mpsfile::String)
The GLPK solver is used to convert the MPS file to a problem with the matrix form
       min  dot(c, x)
subject to bl ≤ A x ≤ bu
           0 ≤ x ≤ xu
"""

function mpstomatrix(mpsfile::String)
    lp_model = MathProgBase.LinearQuadraticModel(ClpSolver())
    # Load the model from .MPS file
    MathProgBase.loadproblem!(lp_model, mpsfile)
    # Transform in vectors
    c = MathProgBase.getobj(lp_model)
    A = MathProgBase.getconstrmatrix(lp_model)
    bl = MathProgBase.getconstrLB(lp_model)
    bu = MathProgBase.getconstrUB(lp_model)
    xl = MathProgBase.getvarLB(lp_model)
    xu = MathProgBase.getvarUB(lp_model)
    return c, A, bl, bu, xl, xu
end

"""
LPtoSTDFormat(c,A,l,u,xlb,xub)
Transforms LP of format
       min  dot(c, x)
subject to l ≦  A x ≦ u
          xl ≤ x ≤ xu
into the standart format
        min  dot(c, x)
subject to  A x = b
          0 ≤ x ≤ xu
"""
function LPtoSTDFormat(c,A,l,u,xlb,xub)
    nrow, ncol = size(A)
    b = Array{Float64,1}()
    # Transform l <= Ax <= u into Ax = b
    delrows = Array{Int64,1}()
    for row = 1:nrow
            if l[row] > u[row]
                throw(error("Problem is infeasible."))
            elseif l[row] == -Inf  && u[row] == Inf #Constraint is always feasible
                push!(delrows,row)
                push!(b,Inf)      #Creates dummy b[row] just to delete at the end.
            elseif l[row] == u[row] # Constraint is l = a'x = u
                    push!(b, l[row])
            elseif l[row] > -Inf && u[row] == Inf #Constraint is  a'x >= l
                ncol += 1
                A = [A spzeros(nrow)]
                A[row,end] = -1.0 # a'x - xs = l
                push!(b, l[row]) # b <- l
                push!(c, 0.) # no cost
                push!(xlb, 0.) #xs >= 0
                push!(xub, Inf) #xs <= Inf
            elseif l[row] == -Inf && u[row] < Inf # Constraint is  a'x <= u
                ncol += 1
                A = [A spzeros(nrow)]
                A[row,end] = 1.0 # a'x + xs = u
                push!(b, u[row]) # b <- u
                push!(c, 0.) # no cost
                push!(xlb, 0.) #xs >= 0
                push!(xub, Inf) #xs <= Inf
            elseif l[row] > -Inf && u[row] < Inf # Constraint is  l <a'x < u.
                A = [A spzeros(nrow)]
                A[row,end] = -1.0 # a'x = xs
                push!(b, 0.) # b <-
                push!(c, 0.)
                push!(xlb, l[row]) # adds variable
                push!(xub, u[row])
                ncol += 1
            end
    end
    A = A[setdiff(1:end,delrows),:]
    b = deleteat!(b,delrows)
    # Transform xlb <= x <= xub  into 0 <= x <= xub
    delcols = Array{Int64,1}()
    for col = 1:ncol
        if xlb[col] > xub[col]
            throw(error("Problem is infeasible."))
        elseif xlb[col] == -Inf  && xub[col] == Inf #Free variable
            A = [A -A[:,col]] #x_i = xp - xm
            push!(c, -c[col]) # adds cost for xm
            xlb[col] = 0.
            push!(xlb, 0.) #xm >= 0
            push!(xub, Inf) #xs <= Inf
        elseif xlb[col] == xub[col] # Constraint is l = x = u
            b -= xlb[col]*A[:,col]
            push!(delcols,col)
        elseif xlb[col] > -Inf && xub[col] <= Inf
            b -= xlb[col]*A[:,col]
            xub[col] -= xlb[col]
            xlb[col] = 0.
        elseif xlb[col] == -Inf && xub[col] < Inf
            c[col] = -c[col]
            A[:,col] = -A[:,col]
            b = b - xub[col]*A[:,col]
            xlb[col] = 0.
            xub[col] = Inf
        end
    end
    A = A[:,setdiff(1:end,delcols)]
    c = c[setdiff(1:end,delcols)]
    xlb = xlb[setdiff(1:end,delcols)]
    xub = xub[setdiff(1:end,delcols)]
    nrow,ncol = size(A)
    return nrow,ncol,c,A,b,xlb,xub
end


"""
FindCircumcentermSet(X)

Finds the Circumcenter of vectors ``x_0,x_1,…,x_m``, columns of matrix ``X``,
as described in [^Behling2018a] and [^Behling2018b].

[^Behling2018]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
Circumcentering the Douglas–Rachford method. Numer. Algorithms. 78(3), 759–776 (2018).
[doi:10.1007/s11075-017-0399-5](https://doi.org/10.1007/s11075-017-0399-5)
[^Behling2018]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
On the linear convergence of the circumcentered-reflection method. Oper. Res. Lett. 46(2), 159-162 (2018).
[doi:10.1016/j.orl.2017.11.018](https://doi.org/10.1016/j.orl.2017.11.018)

"""
    function FindCircumcentermSet(X)
    # Finds the Circumcenter of points X = [X1, X2, X3, ... Xn]
        # println(typeof(X))
        lengthX = length(X)
        if lengthX  == 1
            return X[1]
        elseif lengthX == 2
            return .5*(X[1] + X[2])
        end
        V = []
        b = Float64[]
        # Forms V = [X[2] - X[1] ... X[n]-X[1]]
        # and b = [dot(V[1],V[1]) ... dot(V[n-1],V[n-1])]
        for ind in 2:lengthX
            difXnX1 = X[ind]-X[1]
            push!(V,difXnX1)
            push!(b,dot(difXnX1,difXnX1))
        end

       # Forms Gram Matrix
        dimG = lengthX-1
        G = diagm(b)

        for irow in 1:(dimG-1)
            for icol in  (irow+1):dimG
                G[irow,icol] = dot(V[irow],V[icol])
                G[icol,irow] = G[irow,icol]
            end
        end
        # Can we make this solution faster, or better?
        y = cholesky(G)\b
        CC = X[1]
        for ind in 1:dimG
            CC += .5*y[ind]*V[ind]
        end
        return CC
    end
