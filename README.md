
# DerivativeTests.jl

Since AD is not yet sufficiciently performant for many applications, one
often has to resort to "manual" derivative (gradient, jacobian, hessian, etc)
implementations. This package facilitates testing these implementations using
either AD or FD.

Start using it by
```julia
add DerivativeTests
using DerivativeTests
?fdtest
```

### Example

A test that will pass:
```julia
using DerivativeTests
F(x) = exp(x[1]*x[2]),
dF(x) = [ x[2], x[1] ] * F(x)
fdtest(F, dF, rand(2))
```

A test that will fail:
```julia
using DerivativeTests
F(x) = exp(x[1]*x[2]),
dF(x) = [ x[1], x[2] ] * F(x)
fdtest(F, dF, rand(2))
```
