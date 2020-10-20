function matrix_transform!(B::Array{Complex{Float64},2}, U::AbstractArray{Complex{Float64},2}, A::AbstractArray{Complex{Float64},2})
    # B should be zeroed!
    N=size(A)[1]
    checkbounds(B,1:N,1:N)
    checkbounds(U,1:N,1:N)
    checkbounds(A,1:N,1:N)
    
    # /* assumes that B and A are both Hermitian */

    @inbounds @fastmath for n=1:N
        for q=1:N
            litlsum1 = Complex(0.0)
            for p=1:N
                litlsum1+=U[p,n]*A[q,p]
            end
            for m=1:n
                B[m,n]+=conj(U[q,m])*litlsum1
            end
        end
    end
    for n=1:N
        for m=(n+1):N
            B[m,n]=conj(B[n,m])
        end
    end
end

function matrix_transform!(B::Array{Complex{Float64},2}, U::AbstractArray{Complex{Float64},2}, A::AbstractArray{Complex{Float64},2},u::Integer)
    # B should be zeroed!
    N=size(A)[1]
    checkbounds(B,1:u,1:u)
    checkbounds(U,1:N,1:u)
    checkbounds(A,1:N,1:N)
    
    # /* assumes that B and A are both Hermitian */

    @inbounds @fastmath for n=1:u
        for q=1:N
            litlsum1 = Complex(0.0)
            for p=1:N
                litlsum1+=U[p,n]*A[q,p]
            end
            for m=1:u
                B[m,n]+=conj(U[q,m])*litlsum1
            end
        end
    end
end

export matrix_transform

function matrix_transform(U::AbstractArray{Complex{Float64},2}, A::AbstractArray{Complex{Float64},2})
    B=zeros(Complex{Float64},size(A))
    matrix_transform!(B,U,A)
    B
end

function matrix_transform(U::AbstractArray{Complex{Float64},2}, A::AbstractArray{Complex{Float64},2},u::Integer)
    B=zeros(Complex{Float64},u,u)
    matrix_transform!(B,U,A,u)
    B
end

function matrix_transform3!(B::Array{Complex{Float64},3}, U::AbstractArray{Complex{Float64},2}, Ax::AbstractArray{Complex{Float64},2},Ay::AbstractArray{Complex{Float64},2},Az::AbstractArray{Complex{Float64},2})
    # B should be zeroed!
    N=size(Ax)[1]
    checkbounds(B,1:N,1:N,1:3)
    checkbounds(U,1:N,1:N)
    checkbounds(Ax,1:N,1:N)
    checkbounds(Ay,1:N,1:N)
    checkbounds(Az,1:N,1:N)
    
    # /* assumes that B and A are both Hermitian */

    @inbounds @fastmath for n=1:N
        for q=1:N
            litlsumx = Complex(0.0)
            litlsumy = Complex(0.0)
            litlsumz = Complex(0.0)
            for p=1:N
                litlsumx+=U[p,n]*Ax[q,p]
                litlsumy+=U[p,n]*Ay[q,p]
                litlsumz+=U[p,n]*Az[q,p]
            end
            for m=1:n
                B[m,n,1]+=conj(U[q,m])*litlsumx
                B[m,n,2]+=conj(U[q,m])*litlsumy
                B[m,n,3]+=conj(U[q,m])*litlsumz
            end
        end
    end
    @inbounds for n=1:N
        for m=(n+1):N
            for q=1:3
                B[m,n,q]=conj(B[n,m,q])
            end
        end
    end
end

function matrix_transform3!(B::Array{Complex{Float64},3}, U::AbstractArray{Complex{Float64},2}, Ax::AbstractArray{Complex{Float64},2},Ay::AbstractArray{Complex{Float64},2},Az::AbstractArray{Complex{Float64},2},u::Integer)
    # B should be zeroed!

    checkbounds(B,1:u,1:u,3)
    checkbounds(U,1:14,1:u)
    checkbounds(Ax,1:14,1:14)
    checkbounds(Ay,1:14,1:14)
    checkbounds(Az,1:14,1:14)
    
    # /* assumes that B and A are both Hermitian */

    @inbounds @fastmath for n=1:u
        for q=1:14
            litlsumx = Complex(0.0)
            litlsumy = Complex(0.0)
            litlsumz = Complex(0.0)
            for p=1:14
                litlsumx+=U[p,n]*Ax[q,p]
                litlsumy+=U[p,n]*Ay[q,p]
                litlsumz+=U[p,n]*Az[q,p]
            end
            for m=1:u
                B[m,n,1]+=conj(U[q,m])*litlsumx
                B[m,n,2]+=conj(U[q,m])*litlsumy
                B[m,n,3]+=conj(U[q,m])*litlsumz
            end
        end
    end
end

function matrix_transform3(U::AbstractArray{Complex{Float64},2}, Ax::AbstractArray{Complex{Float64},2},Ay::AbstractArray{Complex{Float64},2},Az::AbstractArray{Complex{Float64},2})
    N=size(Ax)[1]
    B=zeros(Complex{Float64},N,N,3)
    matrix_transform3!(B,U,Ax,Ay,Az)
    B
end

function matrix_transform3(U::AbstractArray{Complex{Float64},2}, Ax::AbstractArray{Complex{Float64},2},Ay::AbstractArray{Complex{Float64},2},Az::AbstractArray{Complex{Float64},2},u::Integer)
    B=zeros(Complex{Float64},u,u,3)
    matrix_transform3!(B,U,Ax,Ay,Az,u)
    B
end
