using LinearAlgebra,DelimitedFiles,StaticArrays
using KDotP
using Plots

pyplot() # or another backend
plot()

kdir=@SVector [0.0,0.0,1.0]

oaxis=0.0:0.005:3.0

a=init_spectrum(oaxis)
a2=init_spectrum2(oaxis)

println("init spectrum")

m=GaAs().parabolic

kx=0.0
ky=0.0

kperp=@SVector [kx,ky,0.0]
println(kperp)
KCMAX=0.5
dkc=1.0/32768
ks=-KCMAX+dkc:dkc:KCMAX
s=calc_c_coeffs(m,kperp,kdir,ks,abstol=1e-6)

l=matrix_element_list(m,kperp,kdir,ks,s)
d=calc_v(l,1:2,1:2)
incr_absorption!(a,m,d)
incr_absorption!(a2,m,d)

#display(plot(a.omega,real(a.v[:,1,1]/2.2918),xlims=(1.3,2.3),label="calc",layout=2))
display(plot(a2.omega,real(a2.v[:,3,3,3,3]),label="calc"))

P0=m.params.P0
μ=0.068
R=3.80998
Eg=m.params.Eg

function one_traj(ħω)
    real(P0^2*μ^0.5/R^0.5/ħω^2*(complex(ħω-Eg))^-0.5)
end

function one_photon_abs(ħω)
    real(P0^2*μ^1.5/R^1.5*2pi/ħω^2*(complex(ħω-Eg))^0.5)
end

function two_photon_one_traj(ħω)
    real(sqrt(R)/(8*μ^0.5*ħω^6)*P0^2*(complex(2ħω-Eg))^0.5)
end

#display(plot!(a.omega,one_traj.(a.omega),label="analytical"))
display(plot!(a2.omega,two_photon_one_traj.(a2.omega),label="analytical"))

if false
for kx=-0.05:0.0005:0.06
    for ky=0.00:0.0005:0.06
        kperp=@SVector [kx,ky,0.0]
        println(kperp)
        KCMAX=0.5
        dkc=1.0/8192
        ks=-KCMAX+dkc:dkc:KCMAX
        s=calc_c_coeffs(m,kperp,kdir,ks,abstol=1e-6)
        if s==nothing
            continue
        end

        l=matrix_element_list(m,kperp,kdir,ks,s)
        d=calc_v(l,KDotP.valence_bands(m),KDotP.conduction_bands(m))
        incr_absorption!(a,m,d)
    end
    #display(plot(a.omega,real(a.v[:,1,1])/2.2918*0.0005^2*2,xlims=(1.3,2.3),label="calc"))
    #display(plot!(a.omega,one_photon_abs.(a.omega)))


end
end
