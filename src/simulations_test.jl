using simulations 

u = simulations.unif(100)
for x in u 
    if x<0 || x>1
        throw(ErrorException("uniform realizations must be between 0 and 1"))
    end
end
e = simulations.inv_exp(u[1], 7)
