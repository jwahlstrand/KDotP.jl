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

function trajectory_intersects_bad(m::Parabolic, kperp, g)
    return false
end
