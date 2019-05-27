module DerivativeTests

using Printf, LinearAlgebra

export fdtest, fdtest_hessian

"""
`fdtest(F, dF, x; kwargs...) -> Bool`

First-order finite-difference test for
   f : ℝⁿ → ℝᵐ
where n, m may be 1, and in this case `Number` rather than `AbstractVector`
arguments and outputs are allows as well.

## Keyword arguments
* `verbose = true` : print test information
* `h0 = 0.1, prange = 2:11` : use finite difference steps `h0.^prange`
* `tol = 1e-3` : required error reduction for succesful test
"""
function fdtest end

# converting scalar inputs into vector inputs
fdtest(F, dF, x::Number; kwargs...) =
   fdtest(y -> F(y[1]), y -> [dF(y[1])], [x]; kwargs...)

# dispatch for different output types
fdtest(F, dF, x; kwargs...) =
   fdtest(F, dF, x, F(x), dF(x); kwargs...)

# convert a scalar output to a vector output
fdtest(F, dF, x::AbstractVector, F0::Number, dF0::AbstractVector; kwargs...) =
   fdtest(x -> [F(x)], Matrix(dF(x)'), x, [F0], Matrix(dF0'); kwargs...)

# the actual finite-difference test implementation 
function fdtest(F, dF, x::AbstractVector, F0::AbstractVector, dF0::AbstractMatrix;
               verbose=true, h0 = 0.1, prange=2:11, tol = 1e-3)
   errors = Float64[]
   dFh = copy(Matrix(dF0))  # the `Matrix` was probably to allow StaticArrays
   @assert size(dFh) == (length(F0), length(x))
   # loop through finite-difference step-lengths
   verbose && @printf("---------|----------- \n")
   verbose && @printf("    h    | error \n")
   verbose && @printf("---------|----------- \n")
   for p = prange
      h = h0^p
      for n = 1:length(x)
         x[n] += h
         dFh[:, n] = (F(x) - F0) / h
         x[n] -= h
      end
      push!(errors, norm(dFh - dF0, Inf))
      verbose && @printf(" %1.1e | %4.2e  \n", h, errors[end])
   end
   verbose && @printf("---------|----------- \n")
   if minimum(errors) <= tol * maximum(errors)
      verbose && @info("passed")
      return true
   end
   @warn("""It seems the finite-difference test has failed, which suggests
      that there might be an inconsistency between the function and gradient
      evaluation. Please double-check this manually / visually. (It is
      also possible that the function being tested is poorly scaled.)""")
   return false
end

end # module
