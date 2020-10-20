export KVector

KDir = SVector{3,Integer}
KVector = SVector{3,Float64}

function three_vector_e(g::KVector)
    if(g[1]==0.0 && g[3]==0.0)
        return KVector(0.0,0.0,1.0)
    end
    e=cross(KVector(0.0,1.0,0.0),g)
    return normalize(e)
end

three_vector_f(g,e) = cross(g,e)

export get_efg,efg_kperp

function get_efg(g)
    e=three_vector_e(g)
    e,three_vector_f(g,e),g
end

efg_kperp(kperp,e,f,g) = kperp[1]*e+kperp[2]*f+kperp[3]*g

function efg_kperp(kperp,kdir)
    e,f,g=get_efg(kdir)
    efg_kperp(kperp,e,f,g)
end

export calc_k

function calc_k(kperp::KVector,kdir::KVector,kc::Float64,pos::Bool)
    kperp + kc*kdir*ifelse(pos,1.0,-1.0)
end

function distance_between_lines(x0,x1,n0,n1)
    if n0 == n1
        return norm(x1-x0)
    end
    n1cn0 = cross(n1, n0)
    return abs(dot(x1-x0,n1cn0))/norm(n1cn0)
end

