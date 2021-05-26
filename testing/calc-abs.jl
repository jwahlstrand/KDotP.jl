using LinearAlgebra,DelimitedFiles,StaticArrays
using KDotP
using Plots

println("done importing stuff")

pyplot() # or another backend

println("pyplot")

plot()

println("plot")

kdir=@SVector [0.0,0.0,1.0]

oaxis=0.0:0.0025:3.0

a=init_spectrum(oaxis)

println("init spectrum")

m=ZnSe().parabolic

KCMAX=0.5
loopRange=0.00:0.0003:0.15

for kx in loopRange
    dkc=1.0/8192
    ks=-KCMAX+dkc:dkc:KCMAX
    for ky in loopRange
        kperp=@SVector [kx,ky,0.0]
        println(kperp)
        s=calc_c_coeffs(m,kperp,kdir,ks,abstol=1e-6)
        if s==nothing
            continue
        end

        l=matrix_element_list(m,kperp,kdir,ks,s)
        d=calc_v(l,KDotP.valence_bands(m),KDotP.conduction_bands(m))
        incr_absorption!(a,m,d)
    end
    display(plot(a.omega,real(a.v[:,1,1]),xlims=(1.3,2.3),label="XX"))
    display(plot!(a.omega,real(a.v[:,2,2]),xlims=(1.3,2.3),label="YY"))
    display(plot!(a.omega,real(a.v[:,3,3]),xlims=(1.3,2.3),label="ZZ"))
end