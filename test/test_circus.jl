include("../src/methods.jl")
#
function test_circus()
    @testset "Simple samples on form min c⋅x  subject to  A x ≦ b" begin
        # Test 1
        c = [-4.; -3]
        A = Matrix([2. 1 2; 3 3 1]')
        b = [4.; 3.; 3.]
        sol = [1.25, .5]
        @show x, f, iter, status = circus(c, A, b)
        @test x ≈ sol
        @test f ≈ dot(sol, c)
        @test status == :Optimal
        # Test 2
        b = [4.; 1; 1]
        A = Matrix([2. 1 2; 3 3 1]')
        c = [-4.; -3]
        sol = [2/5; 1/5]
        @show  x, f, iter, status = circus(c, A, b)
        @test x ≈ sol
        @test f ≈ dot(sol, c)
        @test status == :Optimal
        # Test 3
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
        # Test 4
        c = [0., 0, 1]
        A = [-1. 1 -1; -1 -1 -1; 1 0 0; -1 0. 0.]
        b = [0,0,4.,0.]
        xzero = [1, -1, 0.]
        sol = [4, 0. , -4]
        @show  x, f, iter, status = circus(c, A, b)
        @test x ≈ sol
        @test f ≈ dot(sol, c)
        @test status == :Optimal
    end
end

test_circus()
