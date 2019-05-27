
using DerivativeTests, Test

pass_tests = []

# Test 1: F(x) = exp(x1[1]*x[2])
push!(pass_tests, ( F = x -> exp(x[1]*x[2]),
                   dF = x -> [ x[2], x[1] ] * exp(x[1]*x[2]),
                    x = [1.23, 0.345] ) )

# Test 2: F(x) = ∑ (|x[i]-x[i-1]|^2 + 1)^2
f = r -> 0.25*(r^2+1)^2
df = r -> r^3+r
df2grad = g -> - [g; [0.0] ] + [ [0.0]; g ]
push!(pass_tests, ( F = x -> sum( f(x[i]-x[i-1]) for i = 2:length(x) ),
                   dF = x -> df2grad([ df(x[i]-x[i-1]) for i = 2:length(x) ]),
                    x = rand(10) ) )

# Test 3: F(x) = f(x)
# (this tests the scalar case implementation)
push!(pass_tests, ( F = f, dF = df, x = 0.5678 ))

@info("Running Tests that should pass")
for t in pass_tests
   println(@test fdtest(t.F, t.dF, t.x, verbose=true))
end

# ------------------------------------------------------------

fail_tests = []

# Test 1: F(x) = exp(x1[1]*x[2])
push!(fail_tests, ( F = x -> exp(x[1]*x[2]),
                   dF = x -> [ x[2], x[1] ] * exp(x[2]),
                    x = [1.23, 0.345] ) )

# Test 2: F(x) = ∑ (|x[i]-x[i-1]|^2 + 1)^2
push!(fail_tests, ( F = x -> sum( f(x[i]-x[i-1]) for i = 2:length(x) ),
                   dF = x -> df2grad([ f(x[i]-x[i-1]) for i = 2:length(x) ]),
                    x = rand(10) ) )

# Test 3: F(x) = f(x)
# (this tests the scalar case implementation)
push!(fail_tests, ( F = f, dF = x -> 0.5*df(x), x = 0.5678 ))

@info("Running Tests that should fail")
for t in fail_tests
   println(@test !(fdtest(t.F, t.dF, t.x, verbose=true)))
end
