include("simulations.jl")

n = 1000

u = simulations.unif(n)

for x in u 
    if x<0 || x>1
        throw(ErrorException("uniform realizations must be between 0 and 1"))
    end
end

e = simulations.inv_exp(u, 7)
@assert sum(e.>0) == length(e)


v = n%2==0 ? u : vcat(u, simulations.unif(1)) 
m = length(v)
z = simulations.Box_Muller_method(v[1:Int(m//2)], v[Int(m//2+1):m])[1:n]

g = simulations.Marsaglia_Tsang(1, 3, n=n)
@assert sum(g.>0) == length(g)

b = simulations.beta_from_gamma(2, 89, n=n)
@assert sum(0 .≤b .≤ 1) == length(b)

c = simulations.CRP(n, 10)

dp = simulations.beta_stick_breaking(n, 10)


