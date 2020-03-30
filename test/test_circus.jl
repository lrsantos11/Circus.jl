include("../src/methods.jl")
#
function test_circus()
    @testset "Simple samples on form min c⋅x  subject to  A x ≦ b" begin
        c = [-4.; -3]
        A = Matrix([2. 1 2; 3 3 1]')
        b = [4.; 3.; 3.]
        sol = [1.25, .5]
        @show x, f, iter, status = circus(c, A, b)
        @test x ≈ sol
        @test f ≈ dot(sol, c)
        @test status == :Optimal
        b = [4.; 1; 1]
        A = Matrix([2. 1 2; 3 3 1]')
        c = [-4.; -3]
        sol = [2/5; 1/5]
        @show  x, f, iter, status = circus(c, A, b)
        @test x ≈ sol
        @test f ≈ dot(sol, c)
        @test status == :Optimal
        c = [-3.,-1 , -3]
        A = [2. 1 1; 1 2 3; 2 2 1]
        b = [2., 5, 6]
        sol = [1/5, 0., 8/5]
        A = [A; -I]
        b = [b ;zeros(3)]
        @show  x, f, iter, status = circus(c, A, b)
        @test x ≈ sol
        @test f ≈ dot(sol, c)
        @test status == :Optimal
        c = [4., 5]
        A = [-1. -1; -1 -2; -4 -2; 1 1; 1 -1]
        b = [1., -1., -8, 3, -1]
        sol = [1, 2.]
        A = [A; -I]
        b = [b ;zeros(2)]
        @show  x, f, iter, status = circus(c, A, b)
        @test x ≈ sol
        @test f ≈ dot(sol, c)
        @test status == :Optimal
    end
end

test_circus()
