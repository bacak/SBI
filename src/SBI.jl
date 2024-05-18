"""
Compute a simple blinding index introduced in the paper `A simple blinding index for randomized controlled trials` by Petroff, Bacak, Dagres, Dilk and Wachter, which has been submitted for publication.

Authors: David Petroff and Miroslav Bacak
"""
module SBI

using DataFrames, Distributions, Roots

export blinding_index

default_quantile = quantile(Normal(), 0.975)


"""
    Wilson_CI_z(n_AA::Float64, n_BA::Float64, n_AB::Float64, n_BB::Float64; z::Float64=default_quantile)

Compute the difference of the probabilities p_A-p_B along with the Wilson-CI.

The Wilson-CI is based on Newcombe 1998, PMID: 9595617, Section 2, equations 10.
"""
function Wilson_CI_z(n_AA::Float64, n_BA::Float64, n_AB::Float64, n_BB::Float64; z::Float64=default_quantile)

  m         = n_AA + n_BA
  n         = n_AB + n_BB
  p_A_hat   = n_AA/m
  p_B_hat   = n_AB/n
  theta_hat = p_A_hat - p_B_hat

  # the quadratic equations in Section 2 equations 10 defining l1, l2, u1 and u2 are solved
  l1      = (2*n_AA+z^2)/2/(m+z^2) - z/(m+z^2)*sqrt(z^2/4 + n_AA*(1-p_A_hat))
  l2      = (2*n_AB+z^2)/2/(n+z^2) - z/(n+z^2)*sqrt(z^2/4 + n_AB*(1-p_B_hat))
  u1      = (2*n_AA+z^2)/2/(m+z^2) + z/(m+z^2)*sqrt(z^2/4 + n_AA*(1-p_A_hat))
  u2      = (2*n_AB+z^2)/2/(n+z^2) + z/(n+z^2)*sqrt(z^2/4 + n_AB*(1-p_B_hat))
  delta   = z* sqrt( l1*(1-l1)/(n_AA + n_BA) + u2*(1-u2)/(n_AB + n_BB))
  epsilon = sqrt( (p_A_hat-u1)^2 + (l2-p_B_hat)^2 )
  L       = theta_hat-delta
  U       = theta_hat+epsilon
  
  return DataFrame(est=theta_hat, lwr_ci=L, upr_ci=U)
end

"""
    compute_blinding_index(n_AA::Float64, n_BA::Float64, n_AB::Float64, n_BB::Float64; switch_point::Float64=1E-12, conf_level::Float64=0.95)

Estimate the blinding index, its Wilson CI, p-value dual to the CI and z-value.

See the paper `A simple blinding index for randomized controlled trials` by Petroff, Bacak, Dagres, Dilk and Wachter, which has been submitted for publication.
...
# Arguments
- `n_AA::Float64`: Number of patients in Group A guessing that they are in Group A. A non-negative number, usually an integer.
- `n_BA::Float64`: Number of patients in Group A guessing that they are in Group B. A non-negative number, usually an integer.
- `n_AB::Float64`: Number of patients in Group B guessing that they are in Group A. A non-negative number, usually an integer.
- `n_BB::Float64`: Number of patients in Group B guessing that they are in Group B. A non-negative number, usually an integer.
...
"""
function compute_blinding_index(n_AA::Float64, n_BA::Float64, n_AB::Float64, n_BB::Float64; switch_point::Float64=1E-12, conf_level::Float64=0.95)

  m       = n_AA + n_BA
  n       = n_AB + n_BB
  p_A_hat = n_AA/m
  p_B_hat = n_AB/n
  est     = p_A_hat - p_B_hat    

  n_AA>=0 || throw(DomainError(n_AA, "Argument n_AA cannot be negative."))
  n_BA>=0 || throw(DomainError(n_AA, "Argument n_BA cannot be negative."))
  n_AB>=0 || throw(DomainError(n_AA, "Argument n_AB cannot be negative."))
  n_BB>=0 || throw(DomainError(n_AA, "Argument n_BB cannot be negative."))

  n_AA+n_BA>0 || throw(DomainError([n_AA, n_BA], "n_AA+n_BA must be strictly positive."))
  n_AB+n_BB>0 || throw(DomainError([n_AB, n_BB], "n_AB+n_BB must be strictly positive."))

  switch_point>=0 && switch_point<=1 || throw(DomainError(switch_point, "Argument switch_point must lie in the interval [0, 1]."))
  conf_level>=0   && conf_level<=1   || throw(DomainError(conf_level, "Argument conf_level must lie in the interval [0, 1]."))

  # if the estimate theta_hat is very close to 0 or 1, then the p-value is set by hand
  if abs(est) < switch_point || abs(est) > 1-switch_point
    z_temp  = 0
  end
  if est>=switch_point && est <= 1-switch_point
    function L_zero(z)
      l1_temp = (2*n_AA+z^2)/2/(m+z^2) - z/(m+z^2)*sqrt(z^2/4+n_AA*(1-n_AA/m))
      u2_temp = (2*n_AB+z^2)/2/(n+z^2) + z/(n+z^2)*sqrt(z^2/4+n_AB*(1-n_AB/n))
      L       = l1_temp^2 + u2_temp^2 - 2*n_AA*l1_temp/m - 2*n_AB*u2_temp/n + 2*n_AA*n_AB/m/n
      return L
    end
    z_temp = find_zero(L_zero, (0, 10^3))
  end

  if est <= -switch_point && est >= -1+switch_point
    function U_zero(z)
      l2_temp = (2*n_AA+z^2)/2/(m+z^2) + z/(m+z^2)*sqrt(z^2/4+n_AA*(1-n_AA/m))
      u1_temp = (2*n_AB+z^2)/2/(n+z^2) - z/(n+z^2)*sqrt(z^2/4+n_AB*(1-n_AB/n))
      U       = l2_temp^2 + u1_temp^2 - 2*n_AA*l2_temp/m - 2*n_AB*u1_temp/n + 2*n_AA*n_AB/m/n
      return U 
    end
    z_temp = find_zero(U_zero, (0, 10^3))
  end

  p_value   = 2*(1-cdf.(Normal(), abs(z_temp)))
  wilson_CI = Wilson_CI_z(n_AA, n_BA, n_AB, n_BB, z = quantile(Normal(), 1-(1-conf_level)/2)) 
  df_tmp    = DataFrame(p_value=p_value, z_value=sign(est)*z_temp)
  
  return hcat(wilson_CI, df_tmp)
end 

"""
    blinding_index(freq_table::Matrix{Float64}; switch_point::Float64=1E-12, conf_level::Float64=0.95)

Wrap the function compute_blinding_index.
"""
function blinding_index(freq_table::Matrix{Float64}; switch_point::Float64=1E-12, conf_level::Float64=0.95)
  @assert size(freq_table) == (2, 2)
  n_AA = freq_table[1, 1]
  n_BA = freq_table[2, 1]
  n_AB = freq_table[1, 2]
  n_BB = freq_table[2, 2]
  return compute_blinding_index(n_AA, n_BA, n_AB, n_BB; switch_point, conf_level)
end 


"""
    blinding_index(freq_table::Matrix{Int64}; switch_point::Float64=1E-12, conf_level::Float64=0.95)

Wrap the function compute_blinding_index.
"""
function blinding_index(freq_table::Matrix{Int64}; switch_point::Float64=1E-12, conf_level::Float64=0.95)
  @assert size(freq_table) == (2, 2)
  n_AA = Float64(freq_table[1, 1])
  n_BA = Float64(freq_table[2, 1])
  n_AB = Float64(freq_table[1, 2])
  n_BB = Float64(freq_table[2, 2])
  return compute_blinding_index(n_AA, n_BA, n_AB, n_BB; switch_point, conf_level)
end 


"""
    blinding_index(n_AA::Float64, n_BA::Float64, n_AB::Float64, n_BB::Float64; switch_point::Float64=1E-12, conf_level::Float64=0.95)

Wrap the function compute_blinding_index.
"""
function blinding_index(n_AA::Float64, n_BA::Float64, n_AB::Float64, n_BB::Float64; switch_point::Float64=1E-12, conf_level::Float64=0.95)
  return compute_blinding_index(n_AA, n_BA, n_AB, n_BB; switch_point, conf_level)
end


"""
    blinding_index(n_AA::Integer, n_BA::Integer, n_AB::Integer, n_BB::Integer; switch_point::Float64=1E-12, conf_level::Float64=0.95)

Wrap the function compute_blinding_index.
"""
function blinding_index(n_AA::Integer, n_BA::Integer, n_AB::Integer, n_BB::Integer; switch_point::Float64=1E-12, conf_level::Float64=0.95)
  return compute_blinding_index(Float64(n_AA), Float64(n_BA), Float64(n_AB), Float64(n_BB); switch_point, conf_level)
end

end

