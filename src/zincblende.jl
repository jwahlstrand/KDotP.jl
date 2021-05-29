struct Parameters
    Eg::Float64
    D0::Float64
    E0p::Float64
    D0p::Float64
    db::Float64
    P0::Float64
    Q::Float64
    P0p::Float64
    G1L::Float64
    G2L::Float64
    G3L::Float64
    F::Float64
    Ck::Float64
end

function generate_px(params::Parameters)
    px=zeros(Complex{Float64},14,14)

    px[1,4]=sqrt(3)/sqrt(2)/3*params.P0p
    px[1,13]=-sqrt(3)/3*params.Q
    px[1,14]=sqrt(2)/2*params.Q
    px[2,4]=-sqrt(2)/2*params.P0p
    px[2,12]=-sqrt(3)/3*params.Q
    px[2,14]=sqrt(3)/sqrt(2)/3*params.Q
    px[3,4]=sqrt(3)/3*params.P0p
    px[3,12]=-sqrt(2)/2*params.Q
    px[3,13]=sqrt(3)/sqrt(2)/3*params.Q
    px[4,5]=sqrt(3)/sqrt(2)/3*params.P0
    px[4,6]=-sqrt(2)/2*params.P0
    px[4,7]=sqrt(3)/3*params.P0
    px[5,9]=sqrt(3)/3*params.Q
    px[5,10]=-sqrt(2)/2*params.Q
    px[6,8]=sqrt(3)/3*params.Q
    px[6,10]=-sqrt(3)/sqrt(2)/3*params.Q
    px[7,8]=sqrt(2)/2*params.Q
    px[7,9]=-sqrt(3)/sqrt(2)/3*params.Q
    px[8,11]=sqrt(3)/sqrt(2)/3*params.P0p
    px[9,11]=sqrt(2)/2*params.P0p
    px[10,11]=sqrt(3)/3*params.P0p
    px[11,12]=sqrt(3)/sqrt(2)/3*params.P0
    px[11,13]=sqrt(2)/2*params.P0
    px[11,14]=sqrt(3)/3*params.P0
    return px
end

function generate_py(params::Parameters)
    py=zeros(Complex{Float64},14,14)
    py[1,4]=-sqrt(3)/sqrt(2)/3*params.P0p
    py[1,13]=sqrt(3)/3*params.Q
    py[1,14]=sqrt(2)/2*params.Q
    py[2,4]=-sqrt(2)/2*params.P0p
    py[2,12]=sqrt(3)/3*params.Q
    py[2,14]=-sqrt(3)/sqrt(2)/3*params.Q
    py[3,4]=-sqrt(3)/3*params.P0p
    py[3,12]=-sqrt(2)/2*params.Q
    py[3,13]=-sqrt(3)/sqrt(2)/3*params.Q
    py[4,5]=sqrt(3)/sqrt(2)/3*params.P0
    py[4,6]=sqrt(2)/2*params.P0
    py[4,7]=sqrt(3)/3*params.P0
    py[5,9]=-sqrt(3)/3*params.Q
    py[5,10]=-sqrt(2)/2*params.Q
    py[6,8]=-sqrt(3)/3*params.Q
    py[6,10]=sqrt(3)/sqrt(2)/3*params.Q
    py[7,8]=sqrt(2)/2*params.Q
    py[7,9]=sqrt(3)/sqrt(2)/3*params.Q
    py[8,11]=sqrt(3)/sqrt(2)/3*params.P0p
    py[9,11]=-sqrt(2)/2*params.P0p
    py[10,11]=sqrt(3)/3*params.P0p
    py[11,12]=-sqrt(3)/sqrt(2)/3*params.P0
    py[11,13]=sqrt(2)/2*params.P0
    py[11,14]=-sqrt(3)/3*params.P0

    return py
end

function generate_pz(params::Parameters)
    pz=zeros(Complex{Float64},14,14)

    pz[1,6]=sqrt(3)/3*params.Q
    pz[1,11]=-2*sqrt(3)/sqrt(2)/3*params.P0p
    pz[2,5]=-sqrt(3)/3*params.Q
    pz[2,7]=-2*sqrt(3)/sqrt(2)/3*params.Q
    pz[3,6]=2*sqrt(3)/sqrt(2)/3*params.Q
    pz[3,11]=sqrt(3)/3*params.P0p
    pz[4,8]=2*sqrt(3)/sqrt(2)/3*params.P0p
    pz[4,10]=-sqrt(3)/3*params.P0p
    pz[4,12]=2*sqrt(3)/sqrt(2)/3*params.P0
    pz[4,14]=-sqrt(3)/3*params.P0
    pz[5,11]=-2*sqrt(3)/sqrt(2)/3*params.P0
    pz[7,11]=sqrt(3)/3*params.P0
    pz[8,13]=sqrt(3)/3*params.Q
    pz[9,12]=-sqrt(3)/3*params.Q
    pz[9,14]=-2*sqrt(3)/sqrt(2)/3*params.Q
    pz[10,13]=2*sqrt(3)/sqrt(2)/3*params.Q

    return pz
end

###################
# Trivial 2-band model for checking output with analytical solutions

struct Parabolic <: Model
    params::Parameters
end

nbands(m::Parabolic)=2
valence_bands(m::Parabolic) = 1:1
conduction_bands(m::Parabolic) = 2:2

abstract type Semiconductor <: Model end

struct Semiconductor14 <: Semiconductor
    params::Parameters
    Px::Array{Complex{Float64},2}
    Py::Array{Complex{Float64},2}
    Pz::Array{Complex{Float64},2}
    d2Hdx2::Array{Complex{Float64},2}
    d2Hdy2::Array{Complex{Float64},2}
    d2Hdz2::Array{Complex{Float64},2}
    d2Hdxdy::Array{Complex{Float64},2}
    d2Hdxdz::Array{Complex{Float64},2}
    d2Hdydz::Array{Complex{Float64},2}
    parabolic::Parabolic
end

struct Semiconductor14nr <: Semiconductor
    params::Parameters
    Px::Array{Complex{Float64},2}
    Py::Array{Complex{Float64},2}
    Pz::Array{Complex{Float64},2}
    parabolic::Parabolic
end

function Semiconductor14(params::Parameters)
    px,py,pz=generate_px(params),generate_py(params),generate_pz(params)
    Px,Py,Pz=Hermitian(px),Hermitian(1im*py),Hermitian(pz)
    return Semiconductor14(params,Px,Py,Pz,d2Hdx2(params),d2Hdy2(params),d2Hdz2(params),d2Hdxdy(params),d2Hdxdz(params),d2Hdydz(params),Parabolic(params))
end
Semiconductor14(Eg::AbstractFloat,D0::AbstractFloat,E0p::AbstractFloat,D0p::AbstractFloat,db::AbstractFloat,P0::AbstractFloat,Q::AbstractFloat,P0p::AbstractFloat,G1L::AbstractFloat,G2L::AbstractFloat,G3L::AbstractFloat,F::AbstractFloat,Ck::AbstractFloat)=Semiconductor14(Parameters(Eg,D0,E0p,D0p,db,P0,Q,P0p,G1L,G2L,G3L,F,Ck))

function Semiconductor14nr(params::Parameters)
    px,py,pz=generate_px(params),generate_py(params),generate_pz(params)
    Px,Py,Pz=Hermitian(px),Hermitian(1im*py),Hermitian(pz)
    return Semiconductor14nr(params,Px,Py,Pz,Parabolic(params))
end
Semiconductor14nr(Eg::AbstractFloat,D0::AbstractFloat,E0p::AbstractFloat,D0p::AbstractFloat,db::AbstractFloat,P0::AbstractFloat,Q::AbstractFloat,P0p::AbstractFloat,G1L::AbstractFloat,G2L::AbstractFloat,G3L::AbstractFloat,F::AbstractFloat,Ck::AbstractFloat)=Semiconductor14nr(Parameters(Eg,D0,E0p,D0p,db,P0,Q,P0p,G1L,G2L,G3L,F,Ck))

const ϵ = 0.0

# Hamiltonians and derivatives

function H(m::Semiconductor14nr,k)
    h=zeros(Complex{Float64},14,14)

    E1 = m.params.E0p - m.params.Eg

    G0 = -m.params.Eg-m.params.D0
    G1 = E1+m.params.D0p

    h[1,1]=h[8,8]=G1
    h[2,2]=h[9,9]=G1+ϵ

    h[3,3]=E1
    h[10,10]=E1+ϵ

    h[4,4]=0.0
    h[11,11]=ϵ

    h[5,5]=-m.params.Eg
    h[12,12]=-m.params.Eg+ϵ

    h[6,6]=-m.params.Eg
    h[13,13]=-m.params.Eg+ϵ

    h[7,7]=G0
    h[14,14]=G0+ϵ

    Rk2=R*sum(abs2,k)

    for i=1:14
        h[i,i] += Rk2
    end
    h.=h .+ m.Px .* k[1] .+ m.Py .* k[2] .+ m.Pz .* k[3]

    return Hermitian(h)
end

function dHdx!(h,m::Semiconductor14nr,k)
    copy!(h,m.Px)

    for i=1:14
        h[i,i]+=2*R*k[1]
    end
end

function dHdy!(h,m::Semiconductor14nr,k)
    copy!(h,m.Py)

    for i=1:14
        h[i,i]+=2*R*k[2]
    end
end

function dHdz!(h,m::Semiconductor14nr,k)
    copy!(h,m.Pz)

    for i=1:14
        h[i,i]+=2*R*k[3]
    end
end

nbands(m::Semiconductor14)=14
valence_bands(m::Semiconductor14) =1:6
conduction_bands(m::Semiconductor14) = 7:8

nbands(m::Semiconductor14nr)=14
valence_bands(m::Semiconductor14nr) =1:6
conduction_bands(m::Semiconductor14nr) = 7:8

export H,dHdx,dHdy,dHdz

function H(m::Semiconductor14,k)
    h=zeros(Complex{Float64},14,14)

    Ek=R*sum(abs2,k)
    Ez=R*k[3]^2
    E2zmxy=R*((k[3]+k[1])*(k[3]-k[1])+(k[3]+k[2])*(k[3]-k[2]))

    kx=k[1]
    ky=k[2]
    kz=k[3]

    EP=m.params.P0^2/R
    EQ=m.params.Q^2/R
    g1=m.params.G1L-EP/(3*m.params.Eg)-EQ/(3*m.params.E0p)-EQ/(3*(m.params.E0p+m.params.D0p))
    g2=m.params.G2L-EP/(6*m.params.Eg)+EQ/(6*m.params.E0p)
    g3=m.params.G3L-EP/(6*m.params.Eg)-EQ/(6*m.params.E0p)

    E0 = -m.params.Eg
    E1 = m.params.E0p - m.params.Eg

    G0 = -m.params.Eg-m.params.D0
    G1 = E1+m.params.D0p

    G1p=(G1+E0)/2+sqrt((G1-E0)*(G1-E0)/4-m.params.db*m.params.db/9)
    E0pp=(G1+E0)/2-sqrt((G1-E0)*(G1-E0)/4-m.params.db*m.params.db/9)
    G0p=(G0+E1)/2-sqrt((G0-E1)*(G0-E1)/4-4*m.params.db*m.params.db/9)
    E1p=(G0+E1)/2+sqrt((G0-E1)*(G0-E1)/4-4*m.params.db*m.params.db/9)

    h[1,1]=h[8,8]=G1p+Ek
    h[2,2]=h[9,9]=G1p+Ek+ϵ

    h[3,3]=E1p+Ek
    h[10,10]=E1p+Ek+ϵ

    h[4,4]=2*Ek*m.params.F+Ek
    h[11,11]=2*Ek*m.params.F+Ek+ϵ

    h[5,5]=E0pp-Ek*(g1-g2)-3*Ez*g2
    h[12,12]=E0pp-Ek*(g1-g2)-3*Ez*g2+ϵ

    h[6,6]=E0pp-Ek*(g1+g2)+3*Ez*g2
    h[13,13]=E0pp-Ek*(g1+g2)+3*Ez*g2+ϵ

    h[7,7]=G0p-Ek*g1
    h[14,14]=G0p-Ek*g1+ϵ

    kp = complex(kx,ky)/sqrt(2)
    km=conj(kp)

    h[1,5]=m.params.db/3.0
    h[2,6]=conj(h[1,5])
    h[3,7]=-2*m.params.db/3.0
    h[5,6]=sqrt(3)*R*(kx+ky)*(kx-ky)*g2-m.params.Ck*kz-1im*sqrt(3)*R*kx*ky*2.0*g3
    h[5,7]=sqrt(2)*g2*E2zmxy
    h[5,12]=-sqrt(3)/sqrt(2)*m.params.Ck*kp
    h[5,13]=2*sqrt(3)*sqrt(2)*g3*R*kz*kp-m.params.Ck*km/sqrt(2)
    h[5,14]=6*g3*R*kz*km

    h[6,7]=sqrt(3)*sqrt(2)*R*(kx+ky)*(kx-ky)*g2+1im*sqrt(3)*sqrt(2)*R*kx*ky*2.0*g3
    h[6,12]=2*sqrt(3)*sqrt(2)*g3*R*kz*kp+m.params.Ck*km/sqrt(2)
    h[6,13]=-sqrt(3)/sqrt(2)*m.params.Ck*kp
    h[6,14]=-2*sqrt(3)*g3*R*kz*kp

    h[7,12]=-6*g3*R*kz*km
    h[7,13]=h[6,14]
    h[8,12]=h[1,5]
    h[9,13]=h[1,5]
    h[10,14]=h[2,6]

    h[12,13]=sqrt(3)*R*(kx+ky)*(ky-kx)*g2-m.params.Ck*kz-1im*sqrt(3)*R*kx*ky*2.0*g3
    h[12,14]=h[4,6]
    h[13,14]=sqrt(3)*sqrt(2)*R*(kx+ky)*(ky-kx)*g2+1im*sqrt(3)*sqrt(2)*R*kx*ky*2.0*g3

    h.=h .+ m.Px .* k[1] .+ m.Py .* k[2] .+ m.Pz .* k[3]

    return Hermitian(h)
end


function fill_diags(h,params::Parameters)
    h[1,1]=2*R
    h[2,2]=2*R
    h[8,8]=2*R
    h[9,9]=2*R
    h[3,3]=2*R
    h[10,10]=2*R
    h[4,4]=2*R*2*params.F+2*R
    h[11,11]=h[4,4]
end
fill_diags(h,m::Semiconductor)=fill_diags(h,m.params)

function d2Hdx2(params::Parameters)
    h=zeros(Complex{Float64},14,14)
    EP=params.P0^2/R
    EQ=params.Q^2/R
    g1=params.G1L-EP/(3*params.Eg)-EQ/(3*params.E0p)-EQ/(3*(params.E0p+params.D0p))
    g2=params.G2L-EP/(6*params.Eg)+EQ/(6*params.E0p)

    fill_diags(h,params)
    h[5,5]=-2*R*(g1-g2)
    h[12,12]=h[5,5]
    h[6,6]=-2*R*(g1+g2)
    h[13,13]=h[6,6]
    h[7,7]=-2*R*g1
    h[14,14]=h[7,7]
    h[5,6]=sqrt(3)*2*R*g2
    h[5,7]=-2*R*sqrt(2)*g2
    h[6,7]=sqrt(3)*sqrt(2)*2*R*g2
    h[12,13]=-2*sqrt(3)*R*g2
    h[12,14]=h[5,7]
    h[13,14]=-2*sqrt(3)*sqrt(2)*R*g2
    h
end
d2Hdx2(m::Semiconductor)=d2Hdx2(m.params)

function d2Hdxdy(params::Parameters)
    h=zeros(Complex{Float64},14,14)
    EP=params.P0^2/R
    EQ=params.Q^2/R
    g3=params.G3L-EP/(6*params.Eg)-EQ/(6*params.E0p)

    h[5,6]=-1im*sqrt(3)*R*2*g3
    h[6,7]=1im*sqrt(3)*sqrt(2)*R*2*g3
    h[12,13]=-1im*sqrt(3)*R*2*g3
    h[13,14]=1im*sqrt(3)*sqrt(2)*R*g3*2
    h
end
d2Hdxdy(m::Semiconductor)=d2Hdxdy(m.params)

function d2Hdxdz(params::Parameters)
    h=zeros(Complex{Float64},14,14)
    EP=params.P0^2/R
    EQ=params.Q^2/R
    g3=params.G3L-EP/(6*params.Eg)-EQ/(6*params.E0p)

    h[5,13]=2*sqrt(3)*g3*R
    h[5,14]=6*g3*R/sqrt(2)
    h[6,12]=2*sqrt(3)*g3*R
    h[6,14]=-2*sqrt(3)/sqrt(2)*g3*R
    h[7,12]=-6*g3*R/sqrt(2)
    h[7,13]=h[6,14]
    h
end
d2Hdxdz(m::Semiconductor)=d2Hdxdz(m.params)

function d2Hdy2(params::Parameters)
    h=zeros(Complex{Float64},14,14)
    EP=params.P0^2/R
    EQ=params.Q^2/R
    g1=params.G1L-EP/(3*params.Eg)-EQ/(3*params.E0p)-EQ/(3*(params.E0p+params.D0p))
    g2=params.G2L-EP/(6*params.Eg)+EQ/(6*params.E0p)

    fill_diags(h,params)
    h[5,5]=-2*R*(g1-g2)
    h[12,12]=h[5,5]
    h[6,6]=-2*R*(g1+g2)
    h[13,13]=h[6,6]
    h[7,7]=-2*R*g1
    h[14,14]=h[7,7]

    h[5,6]=-sqrt(3)*2*R*g2
    h[5,7]=-2*R*sqrt(2)*g2
    h[6,7]=-sqrt(3)*sqrt(2)*2*R*g2
    h[12,13]=2*sqrt(3)*R*g2
    h[12,14]=h[5,7]
    h[13,14]=2*sqrt(3)*sqrt(2)*R*g2
    h
end
d2Hdy2(m::Semiconductor)=d2Hdy2(m.params)

function d2Hdydz(params::Parameters)
    h=zeros(Complex{Float64},14,14)
    EP=params.P0^2/R
    EQ=params.Q^2/R
    g3=params.G3L-EP/(6*params.Eg)-EQ/(6*params.E0p)

    h[5,13]=2im*sqrt(3)*g3*R
    h[5,14]=-6im*g3*R/sqrt(2)
    h[6,12]=2im*sqrt(3)*g3*R
    h[6,14]=-2im*sqrt(3)/sqrt(2)*g3*R
    h[7,12]=6im*g3*R/sqrt(2)
    h[7,13]=h[6,14]
    h
end
d2Hdydz(m::Semiconductor)=d2Hdydz(m.params)

function d2Hdz2(params::Parameters)
    h=zeros(Complex{Float64},14,14)
    EP=params.P0^2/R
    EQ=params.Q^2/R
    g1=params.G1L-EP/(3*params.Eg)-EQ/(3*params.E0p)-EQ/(3*(params.E0p+params.D0p))
    g2=params.G2L-EP/(6*params.Eg)+EQ/(6*params.E0p)

    fill_diags(h,params)
    h[5,5]=-2*R*(g1-g2)-6*R*g2
    h[12,12]=h[5,5]
    h[6,6]=-2*R*(g1+g2)+6*R*g2
    h[13,13]=h[6,6]
    h[7,7]=-2*R*g1
    h[14,14]=h[7,7]

    h[5,7]=4*R*sqrt(2)*g2
    h[12,14]=h[5,7]
    h
end
d2Hdz2(m::Semiconductor)=d2Hdz2(m.params)

function dHdx!(h,m::Semiconductor14,k)
    copy!(h,m.Px)

    h[5,12]+=-sqrt(3)*m.params.Ck/2
    h[5,13]+=-m.params.Ck/2
    h[6,12]+=m.params.Ck/2
    h[6,13]+=-sqrt(3)*m.params.Ck/2

    axpy!(k[1],m.d2Hdx2,h)
    axpy!(k[2],m.d2Hdxdy,h)
    axpy!(k[3],m.d2Hdxdz,h)
end

function dHdy!(h,m::Semiconductor14,k)
    copy!(h,m.Py)

    h[5,12]+=-1im*sqrt(3)*m.params.Ck/2
    h[5,13]+=1im*m.params.Ck/2
    h[6,12]+=-1im*m.params.Ck/2
    h[6,13]+=-1im*sqrt(3)*m.params.Ck/2

    axpy!(k[1],m.d2Hdxdy,h)
    axpy!(k[2],m.d2Hdy2,h)
    axpy!(k[3],m.d2Hdydz,h)
end

function dHdz!(h,m::Semiconductor14,k)
    copy!(h,m.Pz)

    h[5,6]+=-m.params.Ck
    h[12,13]+=-m.params.Ck

    axpy!(k[1],m.d2Hdxdz,h)
    axpy!(k[2],m.d2Hdydz,h)
    axpy!(k[3],m.d2Hdz2,h)
end

function H(m::Parabolic,k)
    h=zeros(Complex{Float64},2,2)
    k2=k[1]^2+k[2]^2+k[3]^2
    h[1,1]=-m.params.Eg-R/0.45*k2
    h[2,2]=R/0.08*k2
    Hermitian(h)
end

function dHdx!(h,m::Parabolic,k)
    fill!(h,0.0)

    h[1,2]=m.params.P0
    h[2,1]=m.params.P0

    h[1,1]=-2*R*k[1]/0.45
    h[2,2]=2*R*k[1]/0.08
end

function dHdy!(h,m::Parabolic,k)
    fill!(h,0.0)

    h[1,2]=complex(0.0,m.params.P0)
    h[2,1]=complex(0.0,-m.params.P0)

    h[1,1]=-2*R*k[2]/0.45
    h[2,2]=2*R*k[2]/0.08
end

function dHdz!(h,m::Parabolic,k)
    fill!(h,0.0)

    h[1,2]=m.params.P0
    h[2,1]=m.params.P0

    h[1,1]=-2*R*k[3]/0.45
    h[2,2]=2*R*k[3]/0.08
end

const origin = KVector(0.,0.,0.)

# construct list of lines of degeneracy for both models
degen_list=KVector[]
push!(degen_list,KVector(1.0,0.0,0.0))
push!(degen_list,KVector(0.0,1.0,0.0))
push!(degen_list,KVector(0.0,0.0,1.0))
push!(degen_list,normalize(KVector(1.0,1.0,1.0)))
push!(degen_list,normalize(KVector(-1.0,1.0,1.0)))
push!(degen_list,normalize(KVector(1.0,-1.0,1.0)))
push!(degen_list,normalize(KVector(1.0,1.0,-1.0)))

function trajectory_intersects_bad(m::Semiconductor, kperp, g)
    return any(bad->distance_between_lines(kperp,origin,g,bad)<1e-4,degen_list)
end

function trajectory_intersects_bad(m::Parabolic, kperp, g)
    return false
end

function H(m,k,params::Parameters)
    P0=params.P0
    Q=params.Q
    P0p=params.P0p
    Eg=params.Eg
    E0p=params.E0p
    D0=params.D0
    D0p=params.D0p
    G1L=params.G1L
    G2L=params.G2L
    G3L=params.G3L
    F=params.F
    db=params.db
    Ck=params.Ck
    return H(m,k)
end

#Eg,D0,E0p,D0p,db,P0,Q,P0p,G1L,G2L,G3L,F,Ck
GaAs()=Semiconductor14(1.519,0.341,4.488,0.171,-0.061,10.30,7.70,3.00,7.797,2.458,3.299,-1.055,-3.4)
GaAs_nr()=Semiconductor14nr(1.519,0.341,4.488,0.171,0.0,10.30,7.70,3.00,-1.0,0.0,0.0,0.0,0.0)
ZnSe()=Semiconductor14(2.820,0.403,7.330,0.090,-0.238,10.628,9.845,9.165,4.30,1.14,1.84,0.0,-14.0)
InP()=Semiconductor14(1.424,0.108,4.6,0.50,0.22,8.65,7.24,4.30,5.05,1.6,1.73,0.0,-14)
InSb()=Semiconductor14(0.235,0.803,3.39,0.39,-0.244,9.51,8.22,3.17,40.1,18.1,19.2,0.0,-9.2)
GaSb()=Semiconductor14(0.813,0.75,3.3,0.33,-0.28,9.50,8.12,3.33,13.2,4.4,5.7,0.0,0.43)
