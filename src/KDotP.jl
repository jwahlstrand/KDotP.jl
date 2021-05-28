module KDotP

using LinearAlgebra, OrdinaryDiffEq, StaticArrays

export ufunc, calc_u_coeffs

# Approach is based on Wahlstrand and Sipe, Phys. Rev. B 82, 075206 (2010).
# This paper is referred to in comments below as PRB10.

include("vectors.jl")
include("transforms.jl")

const R = 3.80998 # ħ²/2mₑ in eV ⋅ Å²

export Model

abstract type Model end

function dHdi(m::Model,k,f!)
    n=nbands(m)
    h=zeros(Complex{Float64},n,n)
    f!(h,m,k)
    Hermitian(h)
end

dHdx(m::Model,k)=dHdi(m,k,dHdx!)
dHdy(m::Model,k)=dHdi(m,k,dHdy!)
dHdz(m::Model,k)=dHdi(m,k,dHdz!)

export Zincblende14,Zincblende14nr,Parabolic,Semiconductor,Semiconductor14,Semiconductor14nr,GaAs,ZnSe,GaAs_nr,InP,GaSb,InSb

include("zincblende.jl")

########### Coefficients ODE solving

# We calculate the "coefficients" matrix Dₘₙ(k⟂;kc), defined in PRB10 in Sec.
# IIIB, along a line in the Brillouin zone at a fixed k⟂. Defining g as the
# direction we are calculating along, the line along which matrix elements are
# calculated is k⟂ + kc g. Band energies ħωₘ(k⟂;kc) are calculated
# simultaneously.

# These can be used to find velocity matrix elements at any kc from the velocity
# matrix elements at kc = 0 using Eq. (76). The differential equation we solve
# to find Dₘₙ(k⟂;kc) is Eq. (75).

# In the calculation, we use a real vector u of length 2*N²+N, containing the
# flattened D matrix followed by a vector of band energies.

# The direction along which we calculate is g, and we define perpendicular e and
# f vectors to form a basis. In the calculation, k⟂ is specified in efg basis,
# while g is specified in the xyz basis.

struct ode_params
    m::Model
    kperp::KVector # kperp in the efg basis
    dir::KVector  # g direction
    pos::Bool  # whether we are solving in the +kc or -kc direction
    wc::Matrix{ComplexF64} # pre-allocated matrix
    W::Matrix{ComplexF64}  # pre-allocated matrix for storing the matrix elements in the g direction
    ode_params(m,kperp,dir,pos,n)=new(m,kperp,dir,pos,zeros(Complex{Float64},n,n),zeros(Complex{Float64},n,n))
end

export ode_params

# u holds the coefficients of the matrix elements in the first 2*N^2 elements and the band energies in the last N
function ufunc(du::Vector{Float64},u::Vector{Float64},p::ode_params,kc::Float64)
    k=calc_k(p.kperp,p.dir,kc,p.pos)

    N=size(p.wc)[1]

    # calculate the dH/dq's and generate derivative matrix
    w=p.dir[1]*dHdx(p.m,k)+p.dir[2]*dHdy(p.m,k)+p.dir[3]*dHdz(p.m,k)
    p.wc .= w

    D=reshape(reinterpret(ComplexF64,@view u[1:2*N^2]),(N,N))

    # find Wg at this k
    fill!(p.W,0.0)
    matrix_transform!(p.W,D,p.wc)

    vD=@view du[1:2*N^2] # slicing creates a copy by default

    dD=reshape(reinterpret(ComplexF64,vD),(N,N))
    dD.=0.0

    # calculate change in band energy
    @inbounds for i=1:N
        du[2*N^2+i]=real(p.W[i,i])
    end

    @inbounds for n=1:N
        for q=1:N
            omegadiff = u[2*N^2+n]-u[2*N^2+q]
            if abs(omegadiff)>1e-11
                aa = p.W[q,n]/omegadiff
                for m=1:N
                    dD[m,n]+=aa*D[m,q]
                end
            end
        end
    end

    # all this should be negative if we are going in the -kc direction
    if !p.pos
        @inbounds for i=1:2*N^2+N
            du[i] = -du[i]
        end
    end
end

export initial_u

# The initial matrix Dₘₙ(k⟂;0) and band energies ħωₘ is found by diagonalizing
# the hamiltonian at k⟂.
function initial_u(h)
    n=size(h)[1]
    eigval, eigvec = eigen(h)

    u = zeros(2*n^2+n)

    for i=1:n
        for j=1:n
            u[2*(i-1)+2*n*(j-1)+1]=real(eigvec[i,j])
            u[2*(i-1)+2*n*(j-1)+1+1]=imag(eigvec[i,j])
        end
        u[2*n^2+i]=eigval[i]
    end

    u
end

# maximum kc to calculate - may vary by model so we should allow setting this
const KCMAX=0.5

# this function calculates a unitary matrix that can be used to calculate matrix elements along a direction,
# returning a list of the matrices at values of kc given in ks
function calc_u_coeffs(m::Model, kperp, direction, ks; abstol=5e-7, KCMX=0.5)
    kperp2 = efg_kperp(kperp,direction)

    # if this goes through a degeneracy, just skip it
    if trajectory_intersects_bad(m, kperp2, direction)
        return nothing
    end

    h=H(m,kperp2)

    n=nbands(m)

    u=initial_u(h)

    # make a copy for when we go negative (check if this is really necessary)

    uinit = copy(u)

    pars=ode_params(m,kperp2,direction,true,n)

    # first we solve from 0 to positive values
    prob = ODEProblem(ufunc,u,(0.0,ks[end]),pars)

    dks=ks[2]-ks[1]

    sol_pos = solve(prob,Tsit5(),reltol=5e-7,abstol=abstol,saveat=dks)
    f_pos = [sol_pos[q] for q=1:length(sol_pos)]

    pars=ode_params(m,kperp2,direction,false,n)

    u = copy(uinit)

    # next we solve from 0 to negative values
    prob = ODEProblem(ufunc,u,(0.0,-ks[1]),pars)
    sol_neg = solve(prob,Tsit5(),reltol=5e-7,abstol=abstol,saveat=dks)
    f_neg = [sol_neg[q] for q=2:length(sol_neg)]

    return vcat(reverse(f_neg),f_pos)
end

####### Diagnostics for coefficients ODE solving

export get_energies,get_elements,get_coeffs

"get all band energies from the output of calc_u_coeffs()"
function get_energies(u::AbstractArray,krange)
    [u[q][2*14*14+i] for i=1:14, q=1:length(krange)]
end

export D_from_u

"extract complex D matrix from real u vector"
function D_from_u(u)
    D=reshape(reinterpret(ComplexF64,@view u[1:2*14*14]),(14,14))
end

"extract D row elements from the ith column from the output of calc_u_coeffs()"
function get_elements(u::AbstractArray,krange,i)
    [D_from_u(u[q])[i,j] for j=1:14, q=1:length(krange)]
end

"extract D matrices from the output of calc_u_coeffs()"
function get_coeffs(u::AbstractArray,krange)
    [D_from_u(u[q])[i,j] for i=1:14, j=1:14, q=1:length(krange)]
end

##### Matrix elements

# structure holding band energies and Wx, Wy, and Wz matrix elements at a given
# k = k⟂ + kc g

# W matrix elements are velocity matrix elements -- see Sec. IIIB of PRB10

export matrix_element

struct matrix_element
    k::KVector  # the k vector for this point (in xyz basis)
    kc::Float64 # the kc (k_parallel) value for this point
    ħω::Vector{Float64} # vector of band energies
    W::Array{ComplexF64,3}
end

export matrix_element_from_coeffs

"Calculate Wx, Wy, and Wz matrix elements from a u vector"
function matrix_element_from_coeffs(m::Model,k,u::Vector{Float64},kc::Float64)
    n=nbands(m)
    ħω=u[2*n*n+1:2*n*n+n]

    # 3 allocs here
    D=reshape(reinterpret(ComplexF64,@view u[1:2*n*n]),(n,n))

    W = matrix_transform3(D,dHdx(m,k),dHdy(m,k),dHdz(m,k))

    matrix_element(k,kc,ħω,W)
end

function matrix_element_from_coeffs(m::Model,k,u::Vector{Float64},kc,i::Integer)
    n=nbands(m)
    ħω=u[2*n*n+1:2*n*n+n]

    D=reshape(reinterpret(ComplexF64,@view u[1:2*n*n]),(n,n))

    W = matrix_transform3(D,dHdx(m,k[1]),dHdy(m,k[2]),dHdz(m,k[3]),i)

    matrix_element(k,kc,ħω,W)
end

export matrix_element_list

# get a full list of matrix elements along the line
function matrix_element_list(m::Model,kperp,kdir,kcrange,ca::AbstractArray)
    kperp2 = efg_kperp(kperp,kdir)
    l=Array{matrix_element}(undef,length(kcrange))
    N=nbands(m)
    b1=zeros(ComplexF64,N,N)
    b2=zeros(ComplexF64,N,N)
    b3=zeros(ComplexF64,N,N)
    for q in 1:length(kcrange)
        c=ca[q]
        k=kperp2+kcrange[q]*kdir
        ħω=c[2*N^2+1:2*N^2+N]

        cc=reshape(reinterpret(ComplexF64,@view c[1:2*N^2]),(N,N))
        dHdx!(b1,m,k)
        dHdy!(b2,m,k)
        dHdz!(b3,m,k)

        W = matrix_transform3(cc,b1,b2,b3)

        l[q]=matrix_element(k,kcrange[q],ħω,W)
    end
    l
end

function matrix_element_list(m::Parabolic,kperp,kdir,kcrange,ca::AbstractArray)
    kperp2 = efg_kperp(kperp,kdir)
    l=Array{matrix_element}(undef,length(kcrange))
    N=nbands(m)
    b1=zeros(ComplexF64,N,N)
    b2=zeros(ComplexF64,N,N)
    b3=zeros(ComplexF64,N,N)
    for q in 1:length(kcrange)
        c=ca[q]
        k=kperp2+kcrange[q]*kdir
        ħω=[-1.519-R*norm(k)^2/0.45,R*norm(k)^2/0.08]

        #cc=reshape(reinterpret(Complex{Float64},@view c[1:2*N^2]),(N,N))
        cc=zeros(ComplexF64,N,N)+I
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

        l[q]=matrix_element(k,kcrange[q],ħω,W)
    end
    l
end

const default_Nkc=16384

######### Absorption calculations

# struct containing velocity matrix elements and energy difference between two
# particular bands

struct v_cv
    v::Array{ComplexF64,2} # v matrix elements V_cv (first dim is kc, second dim is [x,y,z])
    ħω::Vector{Float64}          # energy difference ħω_cv for each kc
    dkc::Float64                 # step size in kc, for normalization
end

export calc_v

# Calculates W_{cv}(k⟂;kc)
# For a given ε his can be used to calculate V_{cv}(k⟂;t), and the magnitude of
# this is the same as v_{cv}(k), which is used in Eq. (61) in PRB10 to calculate
# absorption with no DC field

function calc_v(l::Vector{matrix_element},i::Integer,j::Integer;dampfunc=nothing)
    n=length(l)
    v=zeros(ComplexF64,n,3)
    o=zeros(Float64,n)
    q=1
    for me in l
        o[q]=me.ħω[j] - me.ħω[i]
        v[q,:] .= me.W[i,j,:]
        if dampfunc != nothing
            v[q,:] .*= dampfunc(me.kc)
        end
        q=q+1
    end
    dkc=l[2].kc-l[1].kc

    v_cv(v,o,dkc)
end

function calc_v(l::Array{matrix_element,1},ir::UnitRange{Int64},jr::UnitRange{Int64};Nkc=default_Nkc,dampfunc=nothing)
    n=length(l)
    dkc=l[2].kc-l[1].kc
    v=zeros(ComplexF64,n,3,length(ir),length(jr))
    o=zeros(Float64,n,length(ir),length(jr))
    q=1
    for me in l
        if dampfunc != nothing
            damp = dampfunc(me.kc)
        else
            damp = 1.0
        end
        @inbounds for i=ir
            for j=jr
	        o[q,i-ir[1]+1,j-jr[1]+1]=me.ħω[j] - me.ħω[i]

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
            d[(i,j)]=v_cv(v[:,:,i-ir[1]+1,j-jr[1]+1],o[:,i-ir[1]+1,j-jr[1]+1],dkc)
        end
    end
    d
end

####### one-photon absorption

export absorption_spectrum

struct absorption_spectrum
    ħω::Vector{Float64}            # energy in eV
    η::Array{ComplexF64,3}   # transition rate (TODO: specify units)
end

export init_spectrum,incr_absorption!

function init_spectrum(oaxis)
    absorption_spectrum(oaxis,zeros(ComplexF64,length(oaxis),3,3))
end

function Base.:+(a1::absorption_spectrum,a2::absorption_spectrum)
    a=init_spectrum(a1.omega)
    a.η .= a1.η .+ a2.η
    a
end

function scale!(a1::absorption_spectrum,s::Real)
    a1.η .*= s
end

# We calculate the spectrum using a histogram approach. Each wavevector in our
# line contributes to absorption at a particular frequency. We split up
# frequency space into bins and sort the contributions into bins, approximating
# the absorption spectrum.
function incr_absorption!(a::absorption_spectrum,m::Model,d::Dict{Tuple{Int64,Int64},v_cv})
    for vv in valence_bands(m)
        for cc in conduction_bands(m)
            v=d[(vv,cc)]
            # normalizing factor
            # TODO: derive this number from fundamental quantities
            fact=v.dkc/(a.ħω[2]-a.ħω[1])*2.2918
            for q=1:length(v.ħω)
                ħωcv=v.ħω[q]
                ΔE=a.ħω[2]-a.ħω[1]
                # calculate which bin this goes into
                b=Integer(round(Int,ħωcv/ΔE))+1
                if ((b>0) && (b<length(a.ħω)))
                    for j=1:3
                        for i=1:3
                            # Eq. (61) in PRB10
                            γi = v.v[q,i]/ħωcv
                            γj = v.v[q,j]/ħωcv
                            # Multiply by normalizing factor and incremement
                            # this bin in the histogram
                            a.η[b,i,j]+=fact*conj(γi)*γj
                        end
                    end
                end
            end
        end
    end
end

# This just does the matrix element calculation for kperp and kdir and then
# calls incr_absorption!
function incr_absorption!(a::absorption_spectrum,kperp,kdir;abstol=5e-7)
    s=calc_w_phi_coeffs(kperp,kdir,abstol=abstol)
    if s==nothing
        return
    end
    d=calc_v(kperp,kdir,s,1:6,7:8) # valence and conduction bands hard-coded
    incr_absorption!(a,d)
end

# This takes the matrix element list and calculates absorption
function incr_absorption!(a::absorption_spectrum,l::Vector{matrix_element})
    d=calc_v(l,1:6,7:8)  # valence and conduction bands hard-coded
    incr_absorption!(a,d)
end

#### two-photon absorption

# calculate γ, defined in last equation of section 5.1 in notes
function calc_little_gamma2(m::Model,d::Dict{Tuple{Int64,Int64},v_cv},v,c,ωd)
    vcv=d[(v,c)]
    γ=zeros(ComplexF64,length(vcv.ħω),3,3)
    for n=1:nbands(m)
        vcn=d[(n,c)]
        vnv=d[(v,n)]
        for q=1:length(vcv.ħω)
            en=vcv.ħω[q]
            if vcv.ħω[q]>3.0  # no clue why this is here
                continue
            end
            for j=1:3
                for i=1:3
                    γ[q,i,j]+=vcn.v[q,i]*vnv.v[q,j]/(vcn.ħω[q]-vnv.ħω[q]+ωd)+vcn.v[q,j]*vnv.v[q,i]/(vcn.ħω[q]-vnv.ħω[q]-ωd)
                end
            end
        end
    end
    γ .*= (0.5im ./ (vcv.ħω.^2 - ωd^2))  # does not include a whole bunch of other factors
    γ
end

export two_photon_absorption_spectrum, init_spectrum2

struct two_photon_absorption_spectrum
    ħω::Vector{Float64}
    η::Array{ComplexF64,5}
end

function init_spectrum2(oaxis)
    two_photon_absorption_spectrum(oaxis,zeros(ComplexF64,length(oaxis),3,3,3,3))
end

function Base.:+(a1::two_photon_absorption_spectrum,a2::two_photon_absorption_spectrum)
    a=init_spectrum2(a1.ħω)
    a.η .= a1.η .+ a2.η
    a
end

function scale!(a1::two_photon_absorption_spectrum,s::Real)
    a1.η .*= s
end

# calculate the quantity in big parentheses in Eq. (139) in the notes
function incr_absorption!(a::two_photon_absorption_spectrum,m::Model,d::Dict{Tuple{Int64,Int64},v_cv})
    for v=valence_bands(m)
        for c=conduction_bands(m)
            γ=calc_little_gamma2(m,d,v,c,0.0)
            Nkc=size(γ)[1]
            vcv=d[(v,c)]
            fact=4*vcv.dkc/(a.ħω[2]-a.ħω[1])
            for q=1:Nkc
                en=vcv.ħω[q]
                den=a.ħω[2]-a.ħω[1]
                b=Integer(floor(en/2/den))+1
                if ((b>0) && (b<length(a.ħω)))
                    for j=1:3
                        a.η[b,j,j,j,j]+=abs2(γ[q,j,j])*fact
                        for p=1:3
                            if j!=p
                                a.η[b,j,j,p,p]+=conj(γ[q,j,j])*γ[q,p,p]*fact
                                a.η[b,j,p,p,j]+=conj(γ[q,j,p])*γ[q,p,j]*fact
                                a.η[b,j,p,j,p]+=conj(γ[q,j,p])*γ[q,j,p]*fact
                            end
                        end
                    end
                end
            end
        end
    end
end

#######

####### 1+2 photon interference spectrum
export interference_spectrum, init_interference_spectrum

struct interference_spectrum
    ħω::Array{Float64,1}
    η::Array{ComplexF64,4}
end

function init_interference_spectrum(oaxis)
    interference_spectrum(oaxis,zeros(ComplexF64,length(oaxis),3,3,3))
end

function Base.:+(a1::interference_spectrum,a2::interference_spectrum)
    a=init_interference_spectrum(a1.omega)
    a.η .= a1.η .+ a2.η
    a
end

function scale!(a1::interference_spectrum,s::Real)
    a1.η .*= s
end

# calculate the quantity in the last equation of section 10.3.1
function incr_absorption!(a::interference_spectrum,m::Model,d::Dict{Tuple{Int64,Int64},v_cv})
    for v=valence_bands(m)
        for c=conduction_bands(m)
            γ2=calc_little_gamma2(m,d,v,c,0.0)
            Nkc=size(γ)[1]
            vcv=d[(v,c)]
            fact=1*vcv.dkc/(a.ħω[2]-a.ħω[1]) # TODO set this
            for q=1:Nkc
                println(vcv.v[q])
                en=vcv.ħω[q]
                den=a.ħω[2]-a.ħω[1]
                b=Integer(floor(en/2/den))+1
                if ((b>0) && (b<length(a.ħω)))
                    for l=1:3
                        γl = vcv.v[q,i]/en
                        for i=1:3
                            for j=1:3
                                a.η[b,i,j,l]+=conj(γ2[q,i,j])*γl[q,l]*fact #Should this include the complex conjugate
                            end
                        end
                    end
                end
            end
        end
    end
end

######## bundle of lots of spectra

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

# Calculate absorption from one k-space trajectory along direction kdir at
# kperp. This does the full calculation, starting from the Hamiltonian.
function abs_one_traj(m,omega,kperp,kdir)
    a=init_spectrum(omega)
    a2=init_spectrum2(omega)

    KCMAX=0.5
    dkc=1.0/8192
    ks=-KCMAX+dkc:dkc:KCMAX

    s=calc_u_coeffs(m,kperp,kdir,ks,abstol=1e-6)
    if s==nothing
        return spectra(a,a2,0)
    end

    l=matrix_element_list(m,kperp,kdir,ks,s)
    d=calc_v(l,1:14,1:14)
    incr_absorption!(a,m,d)

    incr_absorption!(a2,m,d)

    spectra(a,a2,1)
end

# Sum spectra inside a square box in kperp space.
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
