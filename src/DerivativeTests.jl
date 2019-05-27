module DerivativeTests

using Printf, LinearAlgebra

export fdtest, fdtest_hessian

fdtest(F, dF, x::Number; kwargs...) =
   fdtest(y -> F(y[1]), y -> [dF(y[1])], [x]; kwargs...)


_failure_message() =
   @warn("""It seems the finite-difference test has failed, which suggests
      that there might be an inconsistency between the function and gradient
      evaluation. Please double-check this manually / visually. (It is
      also possible that the function being tested is poorly scaled.)""")


"""
`fdtest(F, dF, x) -> Bool`

First-order finite-difference test for scalar F
"""
function fdtest(F, dF, x::AbstractVector{T};
                verbose=true, h0 = 0.1, prange=2:11
               ) where {T <: AbstractFloat}
   errors = Float64[]
   E = F(x)
   @assert E isa Number  # this version is only valid for scalar-valued F
   dE = dF(x)
   @assert dE isa AbstractVector
   # loop through finite-difference step-lengths
   verbose && @printf("---------|----------- \n")
   verbose && @printf("    h    | error \n")
   verbose && @printf("---------|----------- \n")
   for p in prange
      h = h0^p
      dEh = copy(dE)
      for n = 1:length(dE)
         x[n] += h
         dEh[n] = (F(x) - E) / h
         x[n] -= h
      end
      push!(errors, norm(dE - dEh, Inf))
      @printf(" %1.1e | %4.2e  \n", h, errors[end])
   end
   verbose && @printf("---------|----------- \n")
   if minimum(errors) <= 1e-3 * maximum(errors)
      verbose && @info("passed")
      return true
   else
      _failure_message()
      return false
   end
end


function fdtest_jacobian(F, dF, x::AbstractVector{T};
                         verbose=true, h0 = 0.1, prange=2:11
                        ) where {T <: AbstractFloat}
   errors = Float64[]
   F0 = F(x)
   dF0 = dF(x)
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
   if minimum(errors) <= 1e-3 * maximum(errors)
      verbose && @info("passed")
      return true
   else
      _failure_message()
      return false
   end
end


end # module
