push!(LOAD_PATH,"/home/wahlstrj/Sync/jl/KDotP")
using StaticArrays,BenchmarkTools
using KDotP

kperp=@SVector [0.002,0.004,0.0]
kdir=@SVector [0.0,0.0,1.0]

oaxis=0.0:0.01:3.0
a=init_spectrum(oaxis)

const Nkc=8192
#m=Zincblende14nr()
m=Parabolic()
n=KDotP.nbands(m)

du=zeros(Float64,2*n^2+n)
h=H(m,kperp)
u=initial_c(h)

wc=zeros(Complex{Float64},n,n)
W=zeros(Complex{Float64},n,n)

p=ode_params(m,kperp,kdir,true,n)
@benchmark cfunc(du,u,p,0.0)

const dkc=1.0/Nkc
ks=-0.5+dkc:dkc:0.5
s=calc_c_coeffs(m,kperp,kdir,ks,abstol=1e-6)
s=calc_c_coeffs(m,kperp,kdir,ks,abstol=1e-6)
println("Calculate coeffs:")
@time s=calc_c_coeffs(m,kperp,kdir,ks,abstol=1e-6)

println("Generate matrix element list:")
l=matrix_element_list(m,kperp,kdir,ks,s)
l=matrix_element_list(m,kperp,kdir,ks,s)
@time l=matrix_element_list(m,kperp,kdir,ks,s)

println("Calculate V matrix elements:")
d=calc_v(l,1:14,1:14)
d=calc_v(l,1:14,1:14)
@time d=calc_v(l,1:14,1:14)

println("Calculate one-photon absorption spectra")
incr_absorption!(a,m,d)
incr_absorption!(a,m,d)
@time incr_absorption!(a,m,d)

println("Calculate two-photon absorption spectra")
a2=init_spectrum2(oaxis/2)
incr_absorption!(a2,m,d)
incr_absorption!(a2,m,d)
@time incr_absorption!(a2,m,d)
