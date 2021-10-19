module simulations

using StatsBase

# Linear Congruential Generator with time in nanoseconds as random seed 
function unif(n::Int = 1, m::Int = 2^32, a::Int = 214013, c::Int = 2531011; seed::Integer = 0)
    if seed == 0
        seed = UInt64(time_ns())#*sum(Array{UInt64}(undef, 10))
    else 
        seed = UInt64(seed)
    end
    sample = Float64[]
    for i = 1:n
        seed = (a * seed + c) % m
        push!(sample, seed/m)
    end
    sample        
end

function inv_exp(x::Union{Vector{Float64}, Float64}, lambda::Real = 1)
    return -log.(1 .- x) ./ lambda
end

# sample exponential
function rand_exp(n::Int = 1; lambda::Real = 1)
    u = unif(n)
    return inv_exp(u, lambda)
end

function Box_Muller_method(x::Union{Vector{Float64}, Float64}, y::Union{Vector{Float64}, Float64})
    z_1 = sqrt.(-2 .*log.(x)) .*cos.(2 .*π .*y)    
    z_2 = sqrt.(-2 .*log.(y)) .*cos.(2 .*π .*x)
    return vcat(z_1, z_2)
end

# sample normal
function rand_norm(n::Int = 1; mean::Real = 0, variance::Real = 1)
    m = n%2==0 ? n : n+1 
    u = unif(m)
    z = Box_Muller_method(u[1:Int(m//2)], u[Int(m//2+1):m])
    return variance .* z[1:n] .+ mean
end

# sample gamma 
function Marsaglia_Tsang(n::Int, alpha::Real, beta::Real)
    seed = Int(time_ns())
    d = alpha-1/3
    c = 1/sqrt(9*d)
    sample = []
    while length(sample) != n
        u = unif(3, seed = seed )
        z = Box_Muller_method(u[2], u[3])[1]
        v = (1+c*z)^3
        seed = (seed + Int(time_ns()))%2^32
        if z>-1/c && log(u[1]) < (z^2)/2 + d - d*v + d * log(v)
            push!(sample, d*v)
        end
    end
    return sample ./ beta
end

# sample beta using gamma
function beta_from_gamma(alpha::Real, beta::Real, n::Int =1)
    g_1 = Marsaglia_Tsang(n, alpha, 1)
    g_2 = Marsaglia_Tsang(n, beta, 1)
    b = g_1 ./(g_1 .+ g_2)
    if n==1
        return b[1]
    end
    return b
end

# sample Chinese Restaurant Process
function CRP(n::Int, alpha::Real)
    clients = [1]
    for i=2:n
        weights = ProbabilityWeights(vcat([count(i->(i==x), clients) /(i+alpha-1) for x in unique(clients)], [alpha/(i+alpha-1)]))
        sampler = vcat(unique(clients), [maximum(clients)+1])
        push!(clients, sample(sampler, weights))
    end
    return clients 
end

# sample Dirichlet Process
function beta_stick_breaking(n::Int, alpha::Real)
    probs = [beta_from_gamma(1, alpha)]
    for i = 1:n
        beta = beta_from_gamma(1, alpha)
        push!(probs, beta * cumprod(1 .- probs)[i])
    end
    return probs
end

end # module
