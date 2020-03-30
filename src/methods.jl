export blslp, circumsimplex

include("auxfunctions.jl")

"""
    blslp(c, A, l, u, xl, xu; atol = 1e-5, max_iter = 10000)
Solve the linear program problem on the form
           min  dot(c, x)
           subject to l ≤ A x ≤ u
                      xl ≤ x ≤ xu
where `A` is m × n. `atol` is the tolerance and `max_iter` is the
maximum number of iterations.
"""

function blslp(c, A, l, u, xl, xu; atol = 1e-5, max_iter = 10000)
    b= [u; -l; xu; -xl]
    AA = [A; -A; I; -I]
    return blslp(c, AA, b, atol = atol, max_iter = max_iter)
end

"""
    blslp(c, A, b, xu; atol = 1e-5, max_iter = 10000)
Solve the linear program problem on the form
           min  dot(c, x)
           subject to   A x = b
                      0 ≤ x ≤ xuthrow(ErrorException("test"))
where `A` is m × n. `atol` is the tolerance and `max_iter` is the
maximum number of iterations.
"""
function blslp(c, A, b, xu; atol = 1e-5, max_iter = 10000)
    # Including slacks
    index_slacks =  findall(xu .!=  Inf)
    num_slacks = length(index_slacks)
    num_const, num_var = size(A)
    AA = [A spzeros(num_const,num_slacks);
         sparse(collect(1:num_slacks), index_slacks, ones(Float64,num_slacks),num_slacks,num_var) I]
    bb = [b; xu[index_slacks]]
    cc = [c; zeros(num_slacks)]
    # Computes Solution for the dual
    # min dot(-[b; xu],[y; w])
    # subject to   AA^T[y; w  ≦ [c;0]
    return blslp(-bb, Matrix(AA'), cc, atol = atol, max_iter = max_iter)
end

"""
    blslp(c, A, b; atol = 1e-5, max_iter = 10000)
Solve the linear program problem on the form
           min  dot(c, x)
    subject to  A x ≦ b
where `A` is m × n. `atol` is the tolerance and `max_iter` is the
maximum number of iterations
"""

function blslp(c, A, b; atol = 1e-8, max_iter = 10000)
    # @assert all(b .> zero(T))

    # Set up data structures.
    num_const, num_var = size(A)
    x = zeros(num_var)  # Suppose zero is feasible
    # x[1] = 1
    cnormed = c/norm(c)
    Anormed = copy(A)
    bnormed = copy(b)
    for  i in 1:num_const
        normA = norm(Anormed[i,:])
        Anormed[i,:] /= normA
        bnormed[i] /= normA
    end
    # A = Anormed
    # b = bnormed
    # c = cnormed
    if ~all(A*x.<= b)
        error("x0 is not a feasible starting point")
    end
    # Make Phase I to find a feasible initial solution

    # Using -c as direction to go down
    # d = -c
    # min_ratio = ratiotest(x, A, b, d)
    # Taking the step
    # x += min_ratio*d
    # Begin blslp iterations.
    status = :MaxIter

    iter = 0
    f = dot(c,x)
    while iter <= max_iter
        # Find Circumcenter direction
        d2 = finddirection(x, A, b, Anormed,cnormed)
        # Projects d2 into span{c}
        d1 = d2 - dot(cnormed,d2)*cnormed
        alpha1 = ratiotest(x,A,b,d1)
        alpha2 = ratiotest(x,A,b,d2)
        x1 = x + alpha1*d1
        x2 = x + alpha2*d2
        d3 = finddirection(x2, A, b, Anormed,cnormed)
        # Projects d3 into span{c}
        d3 -= dot(cnormed,d3)*cnormed
        alpha3 = ratiotest(x2,A,b,d3)
        x3 = x2 + alpha3*d3
        # Compute direction
        xm1 = .5*(x + x1)
        xm2 = .5*(x3 + x2)
        d = xm2 - xm1
        alpha = ratiotest(xm1,A,b,d)
        xnew = xm1 + alpha*d
        gaperror = norm(xnew - x)
        fnew = dot(xnew,c)
        # gaperror = norm(fnew - f)
        x = xnew
        iter += 1
        if  gaperror < atol
            # @show x
            status = :Optimal
            # iter
            x = refinesolution(x, A, b, c, num_var,atol)
            f = dot(c,x)
            break
        end
        f = dot(x,c)
    end
    return x,f,iter,status
end

"""
    circumsimplex(c, A, b; atol = 1e-5, max_iter = 10000)
Solve the linear program problem on the form
           min  dot(c, x)
    subject to  A x ≦ b
where `A` is m × n. `atol` is the tolerance and `max_iter` is the
maximum number of iterations
"""

function circumsimplex(c, A, b; xzero = Float64[], atol = 1e-8, max_iter = 10000)
    # @assert all(b .> zero(T))

    # Set up data structures.
    num_const, num_var = size(A)
    isempty(xzero) ?  x = zeros(num_var) : x = xzero
    # Test to verify that xzero is a feasible starting point
    if ~all(A*x - b .<= atol)
        # Make Phase to find a feasible initial solution if necessary
        # by constructing artificial problem
        phaseI = true
        c = [zeros(num_var); 1.]
        maxb = maximum(abs.(b))
        x = [zeros(num_var); maxb]
        b = [b; 0.]
        art_col = spzeros(num_const)
        signb = findall(sign.(b).< 0)
        art_col[signb] = -ones(length(signb))
        A = [A art_col;
            spzeros(num_var)'  -1.]
            num_var +=1
            num_const +=1
    else
        phaseI = false
    end
    cnormed = c/norm(c)
    Anormed = copy(A)
    bnormed = copy(b)
    for  i in 1:num_const
        normA = norm(Anormed[i,:])
        Anormed[i,:] /= normA
        bnormed[i] /= normA
    end
    sortedJ = sortperm((Anormed*cnormed),rev=true)
    # Refine feasible point to find a vertex
    x = refinesolution(x, A, b, c, num_var, atol)
    iter = 0
    f = dot(c,x)
    status = :Optimal

    while iter <= max_iter
        iter += 1
        # Find Circumcenter direction
        d = finddirection(x, A, b, Anormed,cnormed, num_var,sortedJ)
        println(d)
        alpha = ratiotest(x,A,b,d)
        xnew = x + alpha*d
        xnew = refinesolution(xnew, A, b, c, num_var,atol)
        println("Indices das Ativas")
        println(findall(b-A*xnew .<=atol))
        println("Ativas")
        println(Array(xnew))
        gaperror = norm(xnew - x)
        fnew = dot(xnew,c)
        # gaperror = norm(fnew - f)
        x = xnew
        f = dot(c,x)
        if  gaperror < atol
            # @show x
            # status = :Optimal
            # iter
            # x = refinesolution(x, A, b, c, num_var,atol)
            # f = dot(c,x)
            # return x,f,iter,status
            break
        end
    end
    if iter >= max_iter
        status = :MaxIter
    end
    return x,f,iter,status
end
