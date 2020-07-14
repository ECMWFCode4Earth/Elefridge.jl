using Statistics
using Elefridge

i = [0,0,0,0];

round_better(A::Array{Float32}) = abs(mean(A-round(A,7))) <= abs(mean(A-groom(A,7)))

for j in 1:100
    A = rand(Float32,10000000)
    i[1] += round_better(A) ? 1 : 0
    A = A .+ 1
    i[2] += round_better(A) ? 1 : 0
    A = randn(Float32,10000000)
    i[3] += round_better(A) ? 1 : 0
    A = A .+ 1
    i[4] += round_better(A) ? 1 : 0
 end
