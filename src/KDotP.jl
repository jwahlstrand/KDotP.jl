module KDotP

using LinearAlgebra, OrdinaryDiffEq, StaticArrays

export cfunc,calc_w_phi_coeffs, calc_c_coeffs

include("vectors.jl")
include("transforms.jl")

const R = 3.80998 # ħ²/2mₑ in eV ⋅ Å²

export Model

abstract type Model end

function dHdx(m::Model,k)
    n=nbands(m)
    h=zeros(Complex{Float64},n,n)
    dHdx!(h,m,k)
    Hermitian(h)
end

function dHdy(m::Model,k)
    n=nbands(m)
    h=zeros(Complex{Float64},n,n)
    dHdy!(h,m,k)
    Hermitian(h)
end

function dHdz(m::Model,k)
    n=nbands(m)
    h=zeros(Complex{Float64},n,n)
    dHdz!(h,m,k)
    Hermitian(h)
end

export Zincblende14,Zincblende14nr,Parabolic,Semiconductor,Semiconductor14,Semiconductor14nr,GaAs,ZnSe

include("zincblende.jl")
include("parabolic.jl")

########### Coefficients ODE solving

struct ode_params
    m::Model
    kperp::SVector{3,Float64}
    dir::SVector{3,Float64}
    pos::Bool
    wc::Array{Complex{Float64},2}
    W::Array{Complex{Float64},2}
    ode_params(m,kperp,dir,pos,n)=new(m,kperp,dir,pos,zeros(Complex{Float64},n,n),zeros(Complex{Float64},n,n))
end

export ode_params

# u holds the coefficients of the matrix elements in the first 2*N^2 elements and the band energies in the last N
function cfunc(du::Array{Float64,1},u::Array{Float64,1},p::ode_params,kc::Float64)
    k=calc_k(p.kperp,p.dir,kc,p.pos)

    N=size(p.wc)[1]

    # calculate the dH/dq's and generate derivative matrix
    w=p.dir[1]*dHdx(p.m,k)+p.dir[2]*dHdy(p.m,k)+p.dir[3]*dHdz(p.m,k)
    p.wc .= w

    # 3 allocs here
    cc=reshape(reinterpret(Complex{Float64},@view u[1:2*N^2]),(N,N))

    fill!(p.W,0.0)
    matrix_transform!(p.W,cc,p.wc)

    # 1 alloc here
    vcc=@view du[1:2*N^2] # slicing creates a copy by default

    # 2 allocs here
    cdu=reshape(reinterpret(Complex{Float64},vcc),(N,N))
    cdu.=0.0

    @inbounds for i=1:N
        du[2*N^2+i]=real(p.W[i,i])
    end

    @inbounds for n=1:N
        for q=1:N
            omegadiff = u[2*N^2+n]-u[2*N^2+q]
            if abs(omegadiff)>1e-11
                aa = p.W[q,n]/omegadiff
                for m=1:N
                    cdu[m,n]+=aa*cc[m,q]
                end
            end
        end
    end

    if !p.pos
        @inbounds for i=1:2*N^2+N
            du[i] = -du[i]
        end
    end
end

export initial_c

function initial_c(h)
    n=size(h)[1]
    eigval, eigvec = eigen(h)

    c = zeros(2*n^2+n)

    for i=1:n
        for j=1:n
            c[2*(i-1)+2*n*(j-1)+1]=real(eigvec[i,j])
            c[2*(i-1)+2*n*(j-1)+1+1]=imag(eigvec[i,j])
        end
        c[2*n^2+i]=eigval[i]
    end

    c
end

const KCMAX=0.5

# this function calculates a unitary matrix that can be used to calculate matrix elements along a direction,
# returning a function that gives the matrix as a function at a given kc

# confusing but simplifies things:  kperp is in efg basis, direction in xyz basis
# g is the direction of the field.

function calc_w_phi_coeffs(m::Model, kperp, direction; abstol=5e-7, KCMX=0.5)
    kperp2 = efg_kperp(kperp,direction)

    if trajectory_intersects_bad(m, kperp2, direction)
        return nothing
    end

    h=H(m,kperp2)

    n=size(h)[1]

    c=initial_c(h)

    # make a copy for when we go negative

    cinit = copy(c)

    pars=ode_params(kperp2,direction,true,n)

    prob = ODEProblem(cfunc,c,(0.0,KCMX),pars)
    sol_pos = solve(prob,BS3(),reltol=5e-7,abstol=abstol)

    pars=ode_params(kperp2,direction,false,n)

    c = copy(cinit)

    prob = ODEProblem(cfunc,c,(0.0,KCMX),pars)
    sol_neg = solve(prob,BS3(),reltol=5e-7,abstol=abstol)

    function s(kc)
        if kc>=0.0
            sol_pos(kc)
        else
            sol_neg(-kc)
        end
    end

    return s
end

# this function calculates a unitary matrix that can be used to calculate matrix elements along a direction,
# returning a list of the matrices at values of kc given in ks

# confusing but simplifies things:  kperp is in efg basis, direction in xyz basis
# g is the direction of the field.

function calc_c_coeffs(m::Model, kperp, direction, ks; abstol=5e-7, KCMX=0.5)
    kperp2 = efg_kperp(kperp,direction)

    # eventually, split ks into positive and negative parts

    if trajectory_intersects_bad(m, kperp2, direction)
        return nothing
    end

    h=H(m,kperp2)

    n=size(h)[1]

    c=initial_c(h)

    # make a copy for when we go negative

    cinit = copy(c)

    pars=ode_params(m,kperp2,direction,true,n)

    prob = ODEProblem(cfunc,c,(0.0,ks[end]),pars)

    dks=ks[2]-ks[1]

    sol_pos = solve(prob,Tsit5(),reltol=5e-7,abstol=abstol,saveat=dks)
    f_pos = [sol_pos[q] for q=1:length(sol_pos)]

    pars=ode_params(m,kperp2,direction,false,n)

    c = copy(cinit)

    prob = ODEProblem(cfunc,c,(0.0,-ks[1]),pars)
    sol_neg = solve(prob,BS3(),reltol=5e-7,abstol=abstol,saveat=dks)
    f_neg = [sol_neg[q] for q=2:length(sol_neg)]

    return vcat(reverse(f_neg),f_pos)
end

####### Diagnostics for coefficients ODE solving

export get_energies,get_elements,get_coeffs

"get all band energies from the output of calc_c_coeffs()"
function get_energies(c::AbstractArray,krange)
    [c[q][2*14*14+i] for i=1:14, q=1:length(krange)]
end

function get_energies(s::Function,krange)
    [s(k)[2*14*14+i] for i=1:14, k=krange]
end

export C_from_c

"extract complex C matrix from real c vector"
function C_from_c(c)
    cc=reshape(reinterpret(Complex{Float64},@view c[1:2*14*14]),(14,14))
end

"extract C row elements from the ith column from the output of calc_c_coeffs()"
function get_elements(c::AbstractArray,krange,i)
    [C_from_c(c[q])[i,j] for j=1:14, q=1:length(krange)]
end

function get_elements(s::Function,krange,i)
    [C_from_c(s(k))[i,j] for j=1:14, k=krange]
end

"extract C matrices from the output of calc_c_coeffs()"
function get_coeffs(c::AbstractArray,krange)
    [C_from_c(c[q])[i,j] for i=1:14, j=1:14, q=1:length(krange)]
end

function get_coeffs(s::Function,krange)
    [C_from_c(s(k))[i,j] for i=1:14, j=1:14, k=krange]
end

##### Matrix elements

export matrix_element

struct matrix_element
    k::SVector{3,Float64} # the k vector for this point (in efg basis)
    kc::Float64           # the kc (k_parallel) value for this point
    energies::Array{Float64,1}
    W::Array{Complex{Float64},3}
end

export matrix_element_from_coeffs

"Calculate Wx, Wy, and Wz matrix elements from a c vector"
function matrix_element_from_coeffs(m::Model,k,c::Array{Float64,1},kc::Float64)
    energies=c[2*14*14+1:2*14*14+14]

    # 3 allocs here
    cc=reshape(reinterpret(Complex{Float64},@view c[1:2*14*14]),(14,14))

    W = matrix_transform3(cc,dHdx(m,k),dHdy(m,k),dHdz(m,k))

    matrix_element(k,kc,energies,W)
end

function matrix_element_from_coeffs(m::Model,k,c::Array{Float64,1},kc,u::Integer)
    energies=c[2*14*14+1:2*14*14+14]

    cc=reshape(reinterpret(Complex{Float64},@view c[1:2*14*14]),(14,14))

    W = matrix_transform3(cc,dHdx(m,k[1]),dHdy(m,k[2]),dHdz(m,k[3]),u)

    matrix_element(k,kc,energies,W)
end

export matrix_element_list

function matrix_element_list(m,kperp,kdir,kcrange,s::Function)
    kperp2 = efg_kperp(kperp,kdir)
    [matrix_element_from_coeffs(m,kperp2+kc*kdir,s(kc),kc) for kc in kcrange]
end

function matrix_element_list(m::Model,kperp,kdir,kcrange,ca::AbstractArray)
    kperp2 = efg_kperp(kperp,kdir)
    l=Array{matrix_element}(undef,length(kcrange))
    N=nbands(m)
    b1=zeros(Complex{Float64},N,N)
    b2=zeros(Complex{Float64},N,N)
    b3=zeros(Complex{Float64},N,N)
    for q in 1:length(kcrange)
        c=ca[q]
        k=kperp2+kcrange[q]*kdir
        energies=c[2*N^2+1:2*N^2+N]

        cc=reshape(reinterpret(Complex{Float64},@view c[1:2*N^2]),(N,N))
        dHdx!(b1,m,k)
        dHdy!(b2,m,k)
        dHdz!(b3,m,k)

        W = matrix_transform3(cc,b1,b2,b3)

        l[q]=matrix_element(k,kcrange[q],energies,W)
    end
    l
end

function matrix_element_list(m::Parabolic,kperp,kdir,kcrange,ca::AbstractArray)
    kperp2 = efg_kperp(kperp,kdir)
    l=Array{matrix_element}(undef,length(kcrange))
    N=nbands(m)
    b1=zeros(Complex{Float64},N,N)
    b2=zeros(Complex{Float64},N,N)
    b3=zeros(Complex{Float64},N,N)
    for q in 1:length(kcrange)
        c=ca[q]
        k=kperp2+kcrange[q]*kdir
        energies=[-1.519-R*norm(k)^2/0.45,R*norm(k)^2/0.08]

        #cc=reshape(reinterpret(Complex{Float64},@view c[1:2*N^2]),(N,N))
        cc=zeros(Complex{Float64},N,N)+I
        dHdx!(b1,m,kperp2)
        dHdy!(b2,m,kperp2)
        dHdz!(b3,m,kperp2)

        W = matrix_transform3(cc,b1,b2,b3)
        W[1,1,1]=-2*R*k[1]/0.45
        W[2,2,1]=2*R*k[1]/0.08

        W[1,1,2]=-2*R*k[2]/0.45
        W[2,2,2]=2*R*k[2]/0.08

        W[1,1,3]=-2*R*k[3]/0.45
        W[2,2,3]=2*R*k[3]/0.08

        l[q]=matrix_element(k,kcrange[q],energies,W)
    end
    l
end

const default_Nkc=16384

struct v_cv
    v::Array{Complex{Float64},2}
    o::Array{Float64,1}
    max_N::Integer  # maximum N that needs to be kept
    min_N::Integer
    dkc::Float64 # for normalization
end

export calc_v

# calculates v_{cv}(k) [Used in Eq. (61) in PRB (2010)]
function calc_v(kperp,kdir,s::Function,i::Integer,j::Integer;Nkc=default_Nkc)
    dkc=2*KCMAX/Nkc
    denom=KCMAX^4
    x=zeros(Complex{Float64},Nkc)
    y=zeros(Complex{Float64},Nkc)
    z=zeros(Complex{Float64},Nkc)
    o=zeros(Float64,Nkc)
    q=1
    for kc=-KCMAX+dkc:dkc:KCMAX
        if kc>=0.0
            me=matrix_element_from_coeffs(kperp+kc*kdir,s(kc),kc)
        else
            me=matrix_element_from_coeffs(kperp-kc*kdir,s(kc),kc)
        end
        a = me.W[i,j,1]
        b = me.energies[j] - me.energies[i]
        o[q]=b
        damp = exp(-4*abs(kc)^4/denom)

        x[q]=a*damp

        c = me.W[i,j,2]
        y[q]=c*damp

        d = me.W[i,j,3]
        z[q]=d*damp

        q=q+1
    end

    v_cv(x,y,z,o,Nkc/2,0,dkc)
end

function calc_v(l::Array{matrix_element,1},i::Integer,j::Integer)
    n=length(l)
    denom=KCMAX^4
    v=zeros(Complex{Float64},n,3)
    o=zeros(Float64,n)
    q=1
    for me in l
        o[q]=me.energies[j] - me.energies[i]
        v[q,:] .= me.W[i,j,:] .* exp(-4*abs(me.kc)^4/denom)

        q=q+1
    end
    dkc=l[2].kc-l[1].kc

    v_cv(v,o,Nkc/2,0,dkc)
end


function calc_v(kperp,kdir,s::Function,ir::UnitRange{Int64},jr::UnitRange{Int64};Nkc=default_Nkc)
    dkc=2*KCMAX/Nkc
    denom=KCMAX^4
    x=zeros(Complex{Float64},Nkc,length(ir),length(jr))
    y=zeros(Complex{Float64},Nkc,length(ir),length(jr))
    z=zeros(Complex{Float64},Nkc,length(ir),length(jr))
    o=zeros(Float64,Nkc,length(ir),length(jr))
    q=1
    for kc=-KCMAX+dkc:dkc:KCMAX
        if kc>=0.0
            me=matrix_element_from_coeffs(kperp+kc*kdir,s(kc),kc,8)
        else
            me=matrix_element_from_coeffs(kperp-kc*kdir,s(kc),kc,8)
        end
        damp = exp(-4*abs(kc)^4/denom)
        @inbounds for i=ir
            for j=jr
	        o[q,i-ir[1]+1,j-jr[1]+1]=me.energies[j] - me.energies[i]

	        x[q,i-ir[1]+1,j-jr[1]+1]=me.W[i,j,1]*damp
	        y[q,i-ir[1]+1,j-jr[1]+1]=me.W[i,j,2]*damp
	        z[q,i-ir[1]+1,j-jr[1]+1]=me.W[i,j,3]*damp
            end
        end
        q=q+1
    end
    d=Dict{Tuple{Int64,Int64},v_cv}()
    for i=ir
        for j=jr
            d[(i,j)]=v_cv(x[:,i-ir[1]+1,j-jr[1]+1],y[:,i-ir[1]+1,j-jr[1]+1],z[:,i-ir[1]+1,j-jr[1]+1],o[:,i-ir[1]+1,j-jr[1]+1],Nkc/2,0,dkc)
        end
    end
    d
end

function calc_v(l::Array{matrix_element,1},ir::UnitRange{Int64},jr::UnitRange{Int64};Nkc=default_Nkc)
    n=length(l)
    dkc=l[2].kc-l[1].kc
    denom=KCMAX^4
    v=zeros(Complex{Float64},n,3,length(ir),length(jr))
    o=zeros(Float64,n,length(ir),length(jr))
    q=1
    for me in l
        damp = exp(-4*abs(me.kc)^4/denom)
        @inbounds for i=ir
            for j=jr
	        o[q,i-ir[1]+1,j-jr[1]+1]=me.energies[j] - me.energies[i]

                for m=1:3
	            v[q,m,i-ir[1]+1,j-jr[1]+1] = me.W[i,j,m] * damp
                end
            end
        end
        q=q+1
    end
    d=Dict{Tuple{Int64,Int64},v_cv}()
    for i=ir
        for j=jr
            d[(i,j)]=v_cv(v[:,:,i-ir[1]+1,j-jr[1]+1],o[:,i-ir[1]+1,j-jr[1]+1],Nkc/2,0,dkc)
        end
    end
    d
end

####### one-photon absorption

export absorption_spectrum

struct absorption_spectrum
    omega::Array{Float64,1}
    v::Array{Complex{Float64},3}
end

export init_spectrum,incr_absorption!

function init_spectrum(oaxis)
    absorption_spectrum(oaxis,zeros(Complex{Float64},length(oaxis),3,3))
end

function Base.:+(a1::absorption_spectrum,a2::absorption_spectrum)
    a=init_spectrum(a1.omega)
    a.v .= a1.v .+ a2.v
    a
end

function scale!(a1::absorption_spectrum,s::Real)
    a1.v .*= s
end

function incr_absorption!(a::absorption_spectrum,m::Model,d::Dict{Tuple{Int64,Int64},v_cv})
    for vv in valence_bands(m)
        for cc in conduction_bands(m)
            v=d[(vv,cc)]
            fact=v.dkc/(a.omega[2]-a.omega[1])*2.2918
            for q=1:length(v.o)
                en=v.o[q]
                den=a.omega[2]-a.omega[1]
                qq=Integer(round(Int,en/den))+1
                if ((qq>0) && (qq<length(a.omega)))
                    for j=1:3
                        a.v[qq,j,j]+=fact*abs2(v.v[q,j])/en^2
                    end
                end
            end
        end
    end
end

function incr_absorption!(a::absorption_spectrum,kperp,kdir;abstol=5e-7)
    s=calc_w_phi_coeffs(kperp,kdir,abstol=abstol)
    if s==nothing
        return
    end
    d=calc_v(kperp,kdir,s,1:6,7:8)
    incr_absorption!(a,d)
end

function incr_absorption!(a::absorption_spectrum,l::Array{matrix_element,1})
    d=calc_v(l,1:6,7:8)
    incr_absorption!(a,d)
end

#### two-photon absorption

function calc_little_gamma2(m::Model,d::Dict{Tuple{Int64,Int64},v_cv},v,c,omegad)
    vcv=d[(v,c)]
    theta=zeros(Complex{Float64},length(vcv.o),3,3)
    toeV=1240.7/3e5/2/pi
    for n=1:nbands(m)
        vcn=d[(n,c)]
        vnv=d[(v,n)]
        for q=1:length(vcv.o)
            en=vcv.o[q]
            if vcv.o[q]>3.0
                continue
            end
            denom=vcn.o[q]-vnv.o[q]+omegad
            for j=1:3
                theta[q,j,j]+=vcn.v[q,j]*vnv.v[q,j]/denom
                for p=1:3
                    if j!=p
                        theta[q,j,p]+=(vcn.v[q,j]*vnv.v[q,p]+vcn.v[q,p]*vnv.v[q,j])/denom/2
                    end
                end
            end
        end
    end
    theta .*= (0.5im ./ vcv.o.^2)
    theta
end

export two_photon_absorption_spectrum, init_spectrum2

struct two_photon_absorption_spectrum
    omega::Array{Float64,1}
    v::Array{Complex{Float64},5}
end

function init_spectrum2(oaxis)
    two_photon_absorption_spectrum(oaxis,zeros(Complex{Float64},length(oaxis),3,3,3,3))
end

function Base.:+(a1::two_photon_absorption_spectrum,a2::two_photon_absorption_spectrum)
    a=init_spectrum2(a1.omega)
    a.v .= a1.v .+ a2.v
    a
end

function scale!(a1::two_photon_absorption_spectrum,s::Real)
    a1.v .*= s
end

function incr_absorption!(a::two_photon_absorption_spectrum,m::Model,d::Dict{Tuple{Int64,Int64},v_cv})
    for v=valence_bands(m)
        for c=conduction_bands(m)
            theta=calc_little_gamma2(m,d,v,c,0.0)
            Nkc=size(theta)[1]
            vcv=d[(v,c)]
            fact=4*vcv.dkc/(a.omega[2]-a.omega[1])
            for q=1:Nkc
                en=vcv.o[q]
                den=a.omega[2]-a.omega[1]
                qq=Integer(floor(en/2/den))+1
                if ((qq>0) && (qq<length(a.omega)))
                    for j=1:3
                        a.v[qq,j,j,j,j]+=abs2(theta[q,j,j])*fact
                        for p=1:3
                            if j!=p
                                a.v[qq,j,j,p,p]+=conj(theta[q,j,j])*theta[q,p,p]*fact
                                a.v[qq,j,p,p,j]+=conj(theta[q,j,p])*theta[q,p,j]*fact
                                a.v[qq,j,p,j,p]+=conj(theta[q,j,p])*theta[q,j,p]*fact
                            end
                        end
                    end
                end
            end
        end
    end
end

#######

export spectra

struct spectra
    one::absorption_spectrum
    two::two_photon_absorption_spectrum
    n::Integer  # number of spectra used to generate this sum
end

function Base.:+(a1::spectra,a2::spectra)
    new_one=a1.one + a2.one
    new_two=a1.two + a2.two
    spectra(new_one,new_two,a1.n+a2.n)
end

function scale!(a1::spectra,s::Real)
    scale!(a1.one,s)
    scale!(a1.two,s)
end

export abs_one_traj, box_integrate, scale!

function abs_one_traj(m,omega,kperp,kdir)
    a=init_spectrum(omega)
    a2=init_spectrum2(omega)

    KCMAX=0.5
    dkc=1.0/8192
    ks=-KCMAX+dkc:dkc:KCMAX

    s=calc_c_coeffs(m,kperp,kdir,ks,abstol=1e-6)
    if s==nothing
        return spectra(a,a2,0)
    end

    l=matrix_element_list(m,kperp,kdir,ks,s)
    d=calc_v(l,1:14,1:14)
    incr_absorption!(a,m,d)

    incr_absorption!(a2,m,d)

    spectra(a,a2,1)
end

function box_integrate(m,omega,kcent,kwidth,kdir,depth)
    a1=init_spectrum(omega)
    a2=init_spectrum2(omega)
    a=spectra(a1,a2,0)

    println("box: ",kcent)

    kxmin = kcent[1]-kwidth/2
    kxmax = kcent[1]+kwidth/2
    kymin = kcent[2]-kwidth/2
    kymax = kcent[2]+kwidth/2

    if depth==0
        abs_one_traj(m,omega,kcent,kdir)
        scale!(a,kwidth^2)
        return
    end

    n=0
    width = kwidth/depth^2
    for kx = kxmin+width/2:width:kxmax
        for ky = kymin+width/2:width:kymax
            c=SVector(kx, ky, 0)
            #println(c)
            a=a+abs_one_traj(m,omega,c,kdir)
            n=n+1
        end
    end
    scale!(a,width^2)
    if a.n>0
        scale!(a,n/a.n)
    end
    #println(a.n," out of ",n," were used")
    a
end

end # module
