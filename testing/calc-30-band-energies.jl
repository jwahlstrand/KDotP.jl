using LinearAlgebra,DelimitedFiles,StaticArrays
using KDotP
using Plots

pyplot() # or another backend
plot()

m=GaAs30()
ks=-0.02:0.001:0.02

energies=zeros(Float64,length(ks),30)

for i=1:length(ks)
    k=@SVector [ks[i],0.0,0.0] #Just do along the x-direction, although that can easily be changed

    h=H(m,k)

    ens=eigvals(h)
    energies[i,1:end].=ens
end

p=plot(ks,energies)
