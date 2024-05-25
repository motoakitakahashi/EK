# EK

# Julia 1.6.7 
# "Long-term support (LTS) release: v1.6.7 (July 19, 2022)"
# https://julialang.org/downloads/

# install packages
# to use "Pkg.add" I need to declare "using Pkg"
using Pkg

Pkg.add("LinearAlgebra")
using  LinearAlgebra

# define some functions

# number of countries
N = 3

# productivity in countries 
T = [1.0, 2.0, 3.0]

# trade elasticity/shape parameter in Frechet
θ = 4.0

# trade cost matrix 
d = 10.0 * ones(N, N) - Diagonal(fill(9.0, N))
# as in EK 2002, rows are destinations, columns are origins

# labor forces
L = [100.0, 200.0, 300.0]

# the following function gives Φ in EK which is sometimes called the market access term
function market_access(w, d, T, θ)
    output = (d .^ (-θ)) * (T .* (w .^ (-θ)))
    return output
end

# the following function gives π in EK which is called trade shares 
# Btw, in Julia, π and pi is 3.14159...
function trade_shares(w, d, T, Φ, θ, N)
    output = repeat(T' .* (w' .^ (-θ)), N, 1) .* repeat(Φ .^ (-1), 1, N) .* (d .^ (-θ))
    return output
end

# the following function updates wages using labor market clearing and trade balance
function wages(L, π0, w)
    output = ((w .* L)' * π0)' ./ L
    return output
end

tol = 10 ^ (-7)
maxit = 5000

# dampening parameter 
λ = 0.99

# the initial guess for wages 
w_in = [1.0, 1.0, 1.0]

count = 0
dif = 1.0

w = copy(w_in)

while dif > tol && count < maxit 
    Φ = market_access(w, d, T, θ)
    π0 = trade_shares(w, d, T, Φ, θ, N)
    w_new = wages(L, π0, w)
    w_new = w_new ./ fill(w_new[1], N) # the wage in country 1 is the numeraire
    dif = maximum(abs.((w_new - w)./w))
    println("dif: ", dif)
    count = count + 1
    w = λ * w_new + (1-λ) * w
end


w
dif
count