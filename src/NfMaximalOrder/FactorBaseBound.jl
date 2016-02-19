################################################################################
#
#  NfMaximalOrder/FactorbaseBound.jl: Bounds for the factor base for Buchmann
#                                     algorithm
#
# This file is part of Hecke.
#
# Copyright (c) 2015: Claus Fieker, Tommy Hofmann
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#  Copyright (C) 2016 Tommy Hofmann
#
################################################################################

# Computing bounds for the size of the factor base.

################################################################################
#
#  Apres Belabas, Diaz y Diaz, Friedmann
#
################################################################################

# The algorithm is described in
# Belabas, Diaz y Diaz, Friedmann: "Small generators for the ideal class group"

doc"""
***
    factor_base_bound_bdf(O::NfMaximalOrder) -> Int

> Use the algorithm of Belabas, Diaz y Diaz and Friedmann to find $B$ such that
> under GRH the ideal class group of $\mathcal O$ is generated by the prime
> ideals of norm bounded by $B$.
"""
function factor_base_bound_bdf(O::NfMaximalOrder)
  return _factor_base_bound_bdf(O, 100.0, 100.0)
end
  
function _factorbase_bound_bdf_right_side(O::NfMaximalOrder, x0::Float64, D::Dict{Int, Array{Tuple{Int, Int}, 1}})
  K = nf(O) 
  d = degree(K)
  r, s = signature(K)

  curval = Int(ceil(x0))

  R = ArbField(64)

  summ = R(0)

  #curval = 100
  logcurval = log(R(curval))

  # small helper function (is this fast?)
  function comp_summand(p::fmpz, m::Int)
    logp = log(R(p))

    pm2 = R(p)^(R(FlintZZ(m)//FlintZZ(2)))

    secondterm = R(1) - m*logp//logcurval
    
    return logp//pm2 * secondterm
  end

  function comp_summand(p::Int, m::Int)
    return comp_summand(fmpz(p), m)
  end
   
  p = 2

  while p < curval

    max_exp = _max_power_in(p, curval)

    #println("maximal power is $max_exp")

    lP = _prime_dec_cache(O, p, D)

    for P in lP
      Pnorm = fmpz(p)^P[1]
      if Pnorm < curval
        max_exp = _max_power_in(Pnorm, curval)
      
        for m in 1:max_exp
          summand = comp_summand(Pnorm, m)
          summ = summ + summand
        end
      end
    end
    p = next_prime(p)
  end
  
  y = 2*summ - (R(d)*(const_pi(R)^2//R(2)) + r*4*const_catalan(R))//logcurval
  return y::arb
end

function _factor_base_bound_bdf(O::NfMaximalOrder, x0::Float64 = 50.0, ste::Float64 = 20.0)
  K = nf(O) 
  d = degree(K)
  r, s = signature(K)

  R = ArbField(64)

  summ = R(0)

  dec_cache = Dict{Int, Array{Tuple{Int, Int}, 1}}()

  D = log(R(abs(discriminant(O)))) -
      R(degree(K))*(const_euler(R) +
      log(8*const_pi(R))) -
      r * const_pi(R)//R(2)

  y = _factorbase_bound_bdf_right_side(O, x0, dec_cache)

  x1 = 2*x0

  while y < D
    x0 = x1
    x1 = 2*x0
    y = _factorbase_bound_bdf_right_side(O, x1, dec_cache)
  end
      
  dista = abs(x0-x1)

  while !( y > D && dista < ste)
    if y < D 
      x1 = x0 + 3*dista/2
    else
      x1 = x0 - dista/2
    end

    dista = abs(x1-x0)
  
    x0 = x1
    y = _factorbase_bound_bdf_right_side(O, x0, dec_cache)
  end
  @assert ceil(x0) < typemax(Int)
  return Int(ceil(x0))
end

function _prime_dec_cache(O::NfMaximalOrder, p::Int, D::Dict{Int, Array{Tuple{Int, Int}, 1}})
  if haskey(D, p)
    return D[p]
  else
    D[p] = prime_decomposition_type(O, p)
    return D[p]
  end
end

################################################################################
#
#  Apres Bach
#
################################################################################

# The algorithm is described in
# Bach: "Explicit bounds for primality testing and related problems"

doc"""
***
    factor_base_bound_bach(O::NfMaximalOrder) -> Int

> Use the theorem of Bach to find $B$ such that
> under GRH the ideal class group of $\mathcal O$ is generated by the prime
> ideals of norm bounded by $B$.
"""
function factor_base_bound_bach(O::NfMaximalOrder)
  p = 64
  R = ArbField(p)
  if degree(O)==2
    r = ceil(6*log(R(discriminant(O)))^2)
  else
    r = ceil(12*log(R(discriminant(O)))^2)
  end
  (b, x) = unique_integer(r)
  @assert b
  @assert x < typemax(Int)
  return Int(x)::Int
end

################################################################################
#
#  Toplevel function
#
################################################################################

doc"""
***
    factor_base_bound_grh(O::NfMaximalOrder) -> Int

> Returns an integer $B$, such that
> under GRH the ideal class group of $\mathcal O$ is generated by the prime
> ideals of norm bounded by $B$.
"""
function factor_base_bound_grh(O::NfMaximalOrder)
  return min(factor_base_bound_bdf(O), factor_base_bound_bach(O))
end

doc"""
***
    factor_base_bound_minkowski(O::NfMaximalOrder) -> Int

> Returns an integer $B$, such that
> the ideal class group of $\mathcal O$ is generated by the prime
> ideals of norm bounded by $B$.
"""
function factor_base_bound_minkowski(O::NfMaximalOrder)
  r1, r2 = signature(O)
  n = degree(O)
  b = BigInt(round(floor((4/pi)^r2*factorial(BigInt(n))/BigInt(n)^n*sqrt(1.0*abs(discriminant(O))))))
  return fmpz(b)
end

