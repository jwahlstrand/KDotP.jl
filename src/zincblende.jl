
struct Zincblende14nr <: Model end # put parameters in struct?
struct Zincblende14 <: Model
    d2Hdx2::Array{Complex{Float64},2}
    d2Hdy2::Array{Complex{Float64},2}
    d2Hdz2::Array{Complex{Float64},2}
    d2Hdxdy::Array{Complex{Float64},2}
    d2Hdxdz::Array{Complex{Float64},2}
    d2Hdydz::Array{Complex{Float64},2}
end

const P0 = 10.30
const Q = 7.70
const P0p = 3.00
const Eg = 1.519
const E0p = 4.488
const D0 = 0.341
const D0p = 0.171

const G1L = 7.797
const G2L = 2.458
const G3L = 3.299
const F = -1.055
const db = -0.061
const Ck = -0.0034

const epsilon = 0.0

px=zeros(Complex{Float64},14,14)

px[1,4]=sqrt(3)/sqrt(2)/3*P0p
px[1,13]=-sqrt(3)/3*Q
px[1,14]=sqrt(2)/2*Q
px[2,4]=-sqrt(2)/2*P0p
px[2,12]=-sqrt(3)/3*Q
px[2,14]=sqrt(3)/sqrt(2)/3*Q
px[3,4]=sqrt(3)/3*P0p
px[3,12]=-sqrt(2)/2*Q
px[3,13]=sqrt(3)/sqrt(2)/3*Q
px[4,5]=sqrt(3)/sqrt(2)/3*P0
px[4,6]=-sqrt(2)/2*P0
px[4,7]=sqrt(3)/3*P0
px[5,9]=sqrt(3)/3*Q
px[5,10]=-sqrt(2)/2*Q
px[6,8]=sqrt(3)/3*Q
px[6,10]=-sqrt(3)/sqrt(2)/3*Q
px[7,8]=sqrt(2)/2*Q
px[7,9]=-sqrt(3)/sqrt(2)/3*Q
px[8,11]=sqrt(3)/sqrt(2)/3*P0p
px[9,11]=sqrt(2)/2*P0p
px[10,11]=sqrt(3)/3*P0p
px[11,12]=sqrt(3)/sqrt(2)/3*P0
px[11,13]=sqrt(2)/2*P0
px[11,14]=sqrt(3)/3*P0

const Px = Hermitian(px)

py=zeros(Complex{Float64},14,14)

py[1,4]=-sqrt(3)/sqrt(2)/3*P0p
py[1,13]=sqrt(3)/3*Q
py[1,14]=sqrt(2)/2*Q
py[2,4]=-sqrt(2)/2*P0p
py[2,12]=sqrt(3)/3*Q
py[2,14]=-sqrt(3)/sqrt(2)/3*Q
py[3,4]=-sqrt(3)/3*P0p
py[3,12]=-sqrt(2)/2*Q
py[3,13]=-sqrt(3)/sqrt(2)/3*Q
py[4,5]=sqrt(3)/sqrt(2)/3*P0
py[4,6]=sqrt(2)/2*P0
py[4,7]=sqrt(3)/3*P0
py[5,9]=-sqrt(3)/3*Q
py[5,10]=-sqrt(2)/2*Q
py[6,8]=-sqrt(3)/3*Q
py[6,10]=sqrt(3)/sqrt(2)/3*Q
py[7,8]=sqrt(2)/2*Q
py[7,9]=sqrt(3)/sqrt(2)/3*Q
py[8,11]=sqrt(3)/sqrt(2)/3*P0p
py[9,11]=-sqrt(2)/2*P0p
py[10,11]=sqrt(3)/3*P0p
py[11,12]=-sqrt(3)/sqrt(2)/3*P0
py[11,13]=sqrt(2)/2*P0
py[11,14]=-sqrt(3)/3*P0

const Py = Hermitian(1im*py)

pz=zeros(Complex{Float64},14,14)

pz[1,6]=sqrt(3)/3*Q
pz[1,11]=-2*sqrt(3)/sqrt(2)/3*P0p
pz[2,5]=-sqrt(3)/3*Q
pz[2,7]=-2*sqrt(3)/sqrt(2)/3*Q
pz[3,6]=2*sqrt(3)/sqrt(2)/3*Q
pz[3,11]=sqrt(3)/3*P0p
pz[4,8]=2*sqrt(3)/sqrt(2)/3*P0p
pz[4,10]=-sqrt(3)/3*P0p
pz[4,12]=2*sqrt(3)/sqrt(2)/3*P0
pz[4,14]=-sqrt(3)/3*P0
pz[5,11]=-2*sqrt(3)/sqrt(2)/3*P0
pz[7,11]=sqrt(3)/3*P0
pz[8,13]=sqrt(3)/3*Q
pz[9,12]=-sqrt(3)/3*Q
pz[9,14]=-2*sqrt(3)/sqrt(2)/3*Q
pz[10,13]=2*sqrt(3)/sqrt(2)/3*Q

const Pz = Hermitian(pz)

function H(m::Zincblende14nr,k)
    h=zeros(Complex{Float64},14,14)

    E1 = E0p - Eg

    G0 = -Eg-D0
    G1 = E1+D0p

    h[1,1]=h[8,8]=G1
    h[2,2]=h[9,9]=G1+epsilon

    h[3,3]=E1
    h[10,10]=E1+epsilon

    h[4,4]=0.0
    h[11,11]=epsilon

    h[5,5]=-Eg
    h[12,12]=-Eg+epsilon

    h[6,6]=-Eg
    h[13,13]=-Eg+epsilon

    h[7,7]=G0
    h[14,14]=G0+epsilon
    
    Rk2=R*sum(abs2,k)

    for i=1:14
        h[i,i] += Rk2
    end
    h.=h .+ Px .* k[1] .+ Py .* k[2] .+ Pz .* k[3]

    return Hermitian(h)
end

function dHdx!(h,m::Zincblende14nr,k)
    fill!(h,0.0)

    h.+=Px
    for i=1:14
        h[i,i]+=2*R*k[1]
    end
end

function dHdy!(h,m::Zincblende14nr,k)
    fill!(h,0.0)

    h.+=Py
    for i=1:14
        h[i,i]+=2*R*k[2]
    end
end

function dHdz!(h,m::Zincblende14nr,k)
    fill!(h,0.0)

    h.+=Pz
    for i=1:14
        h[i,i]+=2*R*k[3]
    end
end

nbands(m::Zincblende14nr)=14
valence_bands(m::Zincblende14nr) = 1:6
conduction_bands(m::Zincblende14nr) = 7:8

export H,dHdx,dHdy,dHdz

nbands(m::Zincblende14)=14
valence_bands(m::Zincblende14) = 1:6
conduction_bands(m::Zincblende14) = 7:8

function H(m::Zincblende14,k)
    h=zeros(Complex{Float64},14,14)

    Ek=R*sum(abs2,k)
    Ez=R*k[3]^2
    E2zmxy=R*((k[3]+k[1])*(k[3]-k[1])+(k[3]+k[2])*(k[3]-k[2]))

    kx=k[1]
    ky=k[2]
    kz=k[3]
    
    EP=P0^2/R
    EQ=Q^2/R
    g1=G1L-EP/(3*Eg)-EQ/(3*E0p)-EQ/(3*(E0p+D0p))
    g2=G2L-EP/(6*Eg)+EQ/(6*E0p)
    g3=G3L-EP/(6*Eg)-EQ/(6*E0p)

    E0 = -Eg
    E1 = E0p - Eg

    G0 = -Eg-D0
    G1 = E1+D0p

    G1p=(G1+E0)/2+sqrt((G1-E0)*(G1-E0)/4-db*db/9)
    E0pp=(G1+E0)/2-sqrt((G1-E0)*(G1-E0)/4-db*db/9)
    G0p=(G0+E1)/2-sqrt((G0-E1)*(G0-E1)/4-4*db*db/9)
    E1p=(G0+E1)/2+sqrt((G0-E1)*(G0-E1)/4-4*db*db/9)

    h[1,1]=h[8,8]=G1p+Ek
    h[2,2]=h[9,9]=G1p+Ek+epsilon

    h[3,3]=E1p+Ek
    h[10,10]=E1p+Ek+epsilon

    h[4,4]=2*Ek*F+Ek
    h[11,11]=2*Ek*F+Ek+epsilon

    h[5,5]=E0pp-Ek*(g1-g2)-3*Ez*g2
    h[12,12]=E0pp-Ek*(g1-g2)-3*Ez*g2+epsilon

    h[6,6]=E0pp-Ek*(g1+g2)+3*Ez*g2
    h[13,13]=E0pp-Ek*(g1+g2)+3*Ez*g2+epsilon

    h[7,7]=G0p-Ek*g1
    h[14,14]=G0p-Ek*g1+epsilon

    kp = complex(kx,ky)/sqrt(2)
    km=conj(kp)

    h[1,5]=db/3.0
    h[2,6]=conj(h[1,5])
    h[3,7]=-2*db/3.0
    h[5,6]=sqrt(3)*R*(kx+ky)*(kx-ky)*g2-Ck*kz-1im*sqrt(3)*R*kx*ky*2.0*g3
    h[5,7]=sqrt(2)*g2*E2zmxy
    h[5,12]=-sqrt(3)/sqrt(2)*Ck*kp
    h[5,13]=2*sqrt(3)*sqrt(2)*g3*R*kz*kp-Ck*km/sqrt(2)
    h[5,14]=6*g3*R*kz*km

    h[6,7]=sqrt(3)*sqrt(2)*R*(kx+ky)*(kx-ky)*g2+1im*sqrt(3)*sqrt(2)*R*kx*ky*2.0*g3
    h[6,12]=2*sqrt(3)*sqrt(2)*g3*R*kz*kp+Ck*km/sqrt(2)
    h[6,13]=-sqrt(3)/sqrt(2)*Ck*kp
    h[6,14]=-2*sqrt(3)*g3*R*kz*kp

    h[7,12]=-6*g3*R*kz*km
    h[7,13]=h[6,14]
    h[8,12]=h[1,5]
    h[9,13]=h[1,5]
    h[10,14]=h[2,6]
    
    h[12,13]=sqrt(3)*R*(kx+ky)*(ky-kx)*g2-Ck*kz-1im*sqrt(3)*R*kx*ky*2.0*g3
    h[12,14]=h[4,6]
    h[13,14]=sqrt(3)*sqrt(2)*R*(kx+ky)*(ky-kx)*g2+1im*sqrt(3)*sqrt(2)*R*kx*ky*2.0*g3
    
    h.=h .+ Px .* k[1] .+ Py .* k[2] .+ Pz .* k[3]

    return Hermitian(h)
end

function fill_diags(h)
    h[1,1]=2*R
    h[2,2]=2*R
    h[8,8]=2*R
    h[9,9]=2*R
    h[3,3]=2*R
    h[10,10]=2*R
    h[4,4]=2*R*2*F+2*R
    h[11,11]=h[4,4]
end

function d2Hdx2()
    h=zeros(Complex{Float64},14,14)
    EP=P0^2/R
    EQ=Q^2/R
    g1=G1L-EP/(3*Eg)-EQ/(3*E0p)-EQ/(3*(E0p+D0p))
    g2=G2L-EP/(6*Eg)+EQ/(6*E0p)
    
    fill_diags(h)
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

function d2Hdxdy()
    h=zeros(Complex{Float64},14,14)
    EP=P0^2/R
    EQ=Q^2/R
    g3=G3L-EP/(6*Eg)-EQ/(6*E0p)
    
    h[5,6]=-1im*sqrt(3)*R*2*g3
    h[6,7]=1im*sqrt(3)*sqrt(2)*R*2*g3
    h[12,13]=-1im*sqrt(3)*R*2*g3
    h[13,14]=1im*sqrt(3)*sqrt(2)*R*g3*2
    h
end

function d2Hdxdz()
    h=zeros(Complex{Float64},14,14)
    EP=P0^2/R
    EQ=Q^2/R
    g3=G3L-EP/(6*Eg)-EQ/(6*E0p)
    
    h[5,13]=2*sqrt(3)*g3*R
    h[5,14]=6*g3*R/sqrt(2)
    h[6,12]=2*sqrt(3)*g3*R
    h[6,14]=-2*sqrt(3)/sqrt(2)*g3*R
    h[7,12]=-6*g3*R/sqrt(2)
    h[7,13]=h[6,14]
    h
end

function d2Hdy2()
    h=zeros(Complex{Float64},14,14)
    EP=P0^2/R
    EQ=Q^2/R
    g1=G1L-EP/(3*Eg)-EQ/(3*E0p)-EQ/(3*(E0p+D0p))
    g2=G2L-EP/(6*Eg)+EQ/(6*E0p)
    
    fill_diags(h)
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

function d2Hdydz()
    h=zeros(Complex{Float64},14,14)
    EP=P0^2/R
    EQ=Q^2/R
    g3=G3L-EP/(6*Eg)-EQ/(6*E0p)
    
    h[5,13]=2im*sqrt(3)*g3*R
    h[5,14]=-6im*g3*R/sqrt(2)
    h[6,12]=2im*sqrt(3)*g3*R
    h[6,14]=-2im*sqrt(3)/sqrt(2)*g3*R
    h[7,12]=6im*g3*R/sqrt(2)
    h[7,13]=h[6,14]
    h
end

function d2Hdz2()
    h=zeros(Complex{Float64},14,14)
    EP=P0^2/R
    EQ=Q^2/R
    g1=G1L-EP/(3*Eg)-EQ/(3*E0p)-EQ/(3*(E0p+D0p))
    g2=G2L-EP/(6*Eg)+EQ/(6*E0p)
    
    fill_diags(h)
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

function dHdx!(h,m::Zincblende14,k)
    fill!(h,0.0)

    h.+=Px

    h[5,12]+=-sqrt(3)*Ck/2
    h[5,13]+=-Ck/2
    h[6,12]+=Ck/2
    h[6,13]+=-sqrt(3)*Ck/2

    axpy!(k[1],m.d2Hdx2,h)
    axpy!(k[2],m.d2Hdxdy,h)
    axpy!(k[3],m.d2Hdxdz,h)
end

function Zincblende14()
    Zincblende14(d2Hdx2(),d2Hdy2(),d2Hdz2(),d2Hdxdy(),d2Hdxdz(),d2Hdydz())
end

function dHdy!(h,m::Zincblende14,k)
    fill!(h,0.0)

    h.+=Py

    h[5,12]+=-1im*sqrt(3)*Ck/2
    h[5,13]+=1im*Ck/2
    h[6,12]+=-1im*Ck/2
    h[6,13]+=-1im*sqrt(3)*Ck/2

    axpy!(k[1],m.d2Hdxdy,h)
    axpy!(k[2],m.d2Hdy2,h)
    axpy!(k[3],m.d2Hdydz,h)
end

function dHdz!(h,m::Zincblende14,k)
    fill!(h,0.0)

    h.+=Pz

    h[5,6]+=-Ck
    h[12,13]+=-Ck

    axpy!(k[1],m.d2Hdxdz,h)
    axpy!(k[2],m.d2Hdydz,h)
    axpy!(k[3],m.d2Hdz2,h)
end

###################
# Trivial 2-band model for checking output with analytical solutions

struct Parabolic <: Model end

nbands(m::Parabolic)=2
valence_bands(m::Parabolic) = 1:1
conduction_bands(m::Parabolic) = 2:2

function H(m::Parabolic,k)
    h=zeros(Complex{Float64},2,2)
    k2=k[1]^2+k[2]^2+k[3]^2
    h[1,1]=-Eg-R/0.45*k2
    h[2,2]=R/0.08*k2
    Hermitian(h)
end

function dHdx!(h,m::Parabolic,k)
    fill!(h,0.0)

    h[1,2]=P0
    h[2,1]=P0

    h[1,1]=-2*R*k[1]/0.45
    h[2,2]=2*R*k[1]/0.08
end

function dHdy!(h,m::Parabolic,k)
    fill!(h,0.0)

    h[1,2]=complex(0.0,P0)
    h[2,1]=complex(0.0,-P0)

    h[1,1]=-2*R*k[2]/0.45
    h[2,2]=2*R*k[2]/0.08
end

function dHdz!(h,m::Parabolic,k)
    fill!(h,0.0)

    h[1,2]=P0
    h[2,1]=P0

    h[1,1]=-2*R*k[3]/0.45
    h[2,2]=2*R*k[3]/0.08
end

const origin = KVector(0.,0.,0.)
degen_list=KVector[]
push!(degen_list,KVector(1.0,0.0,0.0))
push!(degen_list,KVector(0.0,1.0,0.0))
push!(degen_list,KVector(0.0,0.0,1.0))
push!(degen_list,normalize(KVector(1.0,1.0,1.0)))
push!(degen_list,normalize(KVector(-1.0,1.0,1.0)))
push!(degen_list,normalize(KVector(1.0,-1.0,1.0)))
push!(degen_list,normalize(KVector(1.0,1.0,-1.0)))

function trajectory_intersects_bad(m::Zincblende14, kperp, g)
    return any(bad->distance_between_lines(kperp,origin,g,bad)<1e-4,degen_list)
end

function trajectory_intersects_bad(m::Zincblende14nr, kperp, g)
    return any(bad->distance_between_lines(kperp,origin,g,bad)<1e-4,degen_list)
end

function trajectory_intersects_bad(m::Parabolic, kperp, g)
    return false
end
