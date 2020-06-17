
# Some tests
using Statistics

A = Array{Float32}(1:1e-6:2-1e-6)
Ao = set_one(A,8)
As = shave(A,8)
Ar = round(A,8)
Ag = groom(A,8)

# check bits that are always zero
zero_check = UInt32(0)
for x in As
    global zero_check |= reinterpret(UInt32,x)
end
bitstring(zero_check)

# bias?
mean(As)-mean(A)
mean(Ao)-mean(A)
mean(Ar)-mean(A)
mean(Ag)-mean(A)

# mean error
mean(As-A)
mean(Ao-A)
mean(Ar-A)
mean(Ag-A)

# mean absolute error
mean(abs.(As-A))
mean(abs.(Ao-A))
mean(abs.(Ar-A))
mean(abs.(Ag-A))

# max absolute error
maximum(abs.(As-A))
maximum(abs.(Ao-A))
maximum(abs.(Ar-A))
maximum(abs.(Ag-A))

# max decimal error
maximum(abs.(log10.(As./A)))
maximum(abs.(log10.(Ao./A)))
maximum(abs.(log10.(Ar./A)))
maximum(abs.(log10.(Ag./A)))
