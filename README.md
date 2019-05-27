
# DerivativeTests.jl

Since AD is not yet sufficiciently performant for many applications, one
often has to resort to "manual" derivative (gradient, jacobian, hessian, etc)
implementations. This package facilitates testing these implementations using
either AD or FD.
