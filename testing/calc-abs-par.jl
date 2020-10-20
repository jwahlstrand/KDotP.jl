using Distributed

addprocs(3)

@everywhere using StaticArrays
@everywhere push!(LOAD_PATH,"/home/wahlstrj/repos/KDotP")
@everywhere using KDotP

kdir=@SVector [0.0,0.0,1.0]

m=Zincblende14()

oaxis=0.0:0.01:2.0

const bwidth=0.01
const kpmax=0.1
const depth=3

kperps=[]

for kx=-kpmax+bwidth/2:bwidth:kpmax-bwidth/2
    for ky=bwidth/2:bwidth:kpmax-bwidth
        if sqrt(kx^2+ky^2)<kpmax
            kperp=@SVector [kx,ky,0.0]
            push!(kperps,kperp)
        end
    end
end

#pmap(println,kperps)
#q=pmap(kperp->abs_one_traj(oaxis,kperp,kdir),kperps)

a=@distributed (+) for kperp in kperps
    #abs_one_traj(oaxis,kperp,kdir)
    box_integrate(m,oaxis,kperp,bwidth,kdir,depth)
end

scale!(a,2.0) # to account for the fact that we integrated over half the BZ

println("all done")
