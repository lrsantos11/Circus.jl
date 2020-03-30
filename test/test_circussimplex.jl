include("../src/methods.jl")

function test_circussimplex()
    # Needs to make MPS available
    # @testset "Test with MPS" begin
    #     mpsfile = "mps/afiro.mps"
    #     num_const, num_var,cc, AA, bb, xu = mpstomatrix(mpsfile)
    #     # Writing the Dual Problem
    #     J = findall(xu .!= Inf)
    #     num_bounds = length(J)
    #     IJ = spzeros(num_bounds,num_var)
    #     for i = 1:num_bounds
    #         IJ[i,J[i]] = 1.
    #     end

    #     A  = [AA' IJ';
    #           spzeros(num_bounds,num_const) I]
    #     c = -[bb; xu[J]]
    #     b = [cc; zeros(num_bounds)]


    #     model = Model(with_optimizer(GLPK.Optimizer,method=:InteriorPoint,msg_lev=GLPK.MSG_ON,))
    #     @variable(model, x[1:(num_const+num_bounds)])
    #     @objective(model, Min, dot(zeros(num_const+num_bounds),x))
    #     @constraint(model,A*x .<= b)
    #     optimize!(model)
    #     println(JuMP.termination_status(model))
    #     println(JuMP.primal_status(model))
    #     println(JuMP.objective_value(model))
    #     xzero = value.(x)
    #     x, f, iter, status = circussimplex(c, A, b, xzero = xzero,  max_iter = 100)

    #     # @test x ≈ sol
    #     # @test f ≈ dot(sol, c)
    #     @test status == :Optimal
    # end
end

test_circussimplex()