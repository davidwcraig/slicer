# import Tmp # import module of slicerm
# commented out to use from REPL

using LinearAlgebra # to use dot fun

function g2(x, xa)
    return -0.5*dot(x,x)
end

res = Tmp.slicerm(g2, 2, [0.0 0.1], 0, N=500) # run a 500 sample MCMC



