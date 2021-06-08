using LinearAlgebra
struct Ket
    factor
    type
    spin::Bool #true if up
end

#Holds an array of kets in which each element is a term
struct Terms
    terms::Array{Ket,1}
end
#Initialize a Terms with only one ket
Terms(ket::Ket)=Terms([ket])
#Used to create a zeros(Terms) matrix
Base.zero(Terms)=Terms([])

struct Spin
    spin::Bool
end
sup=Spin(true)
sdown=Spin(false)

#Define an array of functions which generate orthogonal states, each labelled with a name, which defaults to an integer count
genBasis(names::Array=[])=[function (spin::Bool=true) Ket(1,names[i],spin) end for i=1:length(names)]
genBasis(numStates::Integer)=genBasis(Array(1:numStates))

#For printing purposes
Base.string(ket::Ket)=string((typeof(ket.factor)== Integer && ket.factor==1 ? "" : ket.factor), "| ", ket.type,(ket.spin ? "↑" : "↓"), " ⟩")
Base.string(terms::Terms)=begin
    if length(terms.terms)==0
        return ""
    end
    str=string(terms.terms[1])
    for term in terms.terms[2:end]
        str=string(str, " + ",term)
    end
    str
end
Base.show(io::IO,ket::Ket)=print(io,string(ket))
Base.show(io::IO,terms::Terms)=print(io,string(terms))

#Functions to create a ket
x,y,z,X,Y,Z,s=genBasis(['x','y','z','X','Y','Z','s'])

##Group algebra rules

# a::Ket * b::Ket is equivalent to ⟨a|b⟩
Base.:*(a::Ket,b::Ket)=a.type==b.type && a.spin==b.spin ? conj(a.factor)*b.factor : 0
#When multiplying a ket by a factor, multiply the ket's factor field by that factor
Base.:*(ket::Ket,fact)=Ket(ket.factor*fact,ket.type,ket.spin)
Base.:*(fact,ket::Ket)=ket*fact
#Adding kets returns a new ket if they're the same ket and the same spin, otherwise returns a Terms
Base.:+(a::Ket,b::Ket)=begin
    if a.type==b.type && a.spin==b.spin
        return Ket(a.factor+b.factor,a.type,a.spin)
    else
        return Terms([a,b])
    end
end
#When adding a Ket to a Terms, either append a new term or add the Ket's factor to an existing term
Base.:+(a::Ket,b::Terms)=begin
    for i in 1:length(b.terms)
        if b.terms[i].type==a.type && b.terms[i].spin==a.spin
            new_terms=copy(b.terms)
            new_terms[i]=b.terms[i]+a
            return Terms(new_terms)
        end
    end
    return Terms([b.terms;a])
end
Base.:+(b::Terms,a::Ket)=a+b
#Adding Terms
Base.:+(a::Terms,b::Terms)=begin
    sum=Terms(a.terms)
    for ket in b.terms
        sum+=ket
    end
    sum
end
#Multiplying Terms by a factor is just distribution
Base.:*(fact,terms::Terms)=Terms(fact .* terms.terms)
Base.:*(terms::Terms,fact)=fact*terms
#Multiplying Terms by a Ket applies the bra to each term and sums up
Base.:*(fact::Ket,terms::Terms)=begin
    sum=0
    for ket in terms.terms
        sum+=fact*ket
    end
    sum
end
Base.:*(terms::Terms,fact::Ket)=fact*Terms
#Multiplying two Terms means distribution
Base.:*(a::Terms,b::Terms)=begin
    new_terms=0
    for term in a.terms
        new_terms+=term*b
    end
    new_terms
end

Base.sum(terms::Terms)=begin
    sum=zero(Terms)
    for term in terms.terms
        sum+=term
    end
    sum
end
Base.sum(terms::Array{Ket,1})=begin
    sum=zero(Terms)
    for term in terms
        sum+=term
    end
    sum
end
Base.:*(a::Ket,spin::Spin)=Ket(a.factor,a.type,spin.spin)
Base.:*(spin::Spin,a::Ket)=a*spin
Base.:*(a::Terms,spin::Spin)=Terms([term*spin for term in a.terms])
Base.:*(spin::Spin,a::Terms)=a*spin
Base.:*(a::Ket,spin::Bool)=Ket(a.factor,a.type,spin)
Base.:*(spin::Bool,a::Ket)=a*spin
Base.:*(a::Terms,spin::Bool)=Terms([term*spin for term in a.terms])
Base.:*(spin::Bool,a::Terms)=a*spin

#Define the basis states in the order given by the Mathematica notebook
basis=zeros(Terms,14)
basis[1]=-sqrt(2/3)*z()+(1/sqrt(6))*(x(false)+1im*y(false))
basis[2]=-sqrt(1/2)*(x(false)+-1im*y(false))
basis[3]=sqrt(1/3)*z()+sqrt(1/3)*(x(false)+1im*y(false))
basis[4]=Terms([1im*s(false)])
basis[5]=-sqrt(2/3)*Z()+sqrt(1/6)*(X(false)+1im*Y(false))
basis[6]=-sqrt(1/2)*(X(false)+-1im*Y(false))
basis[7]=sqrt(1/3)*Z()+sqrt(1/3)*(X(false)+1im*Y(false))
basis[8]=sqrt(2/3)*z(false)+sqrt(1/6)*(x()+-1im*y())
basis[9]=sqrt(1/2)*(x()+1im*y())
basis[10]=-sqrt(1/3)*z(false)+sqrt(1/3)*(x()+-1im*y())
basis[11]=Terms([1im*s()])
basis[12]=sqrt(2/3)*Z(false)+sqrt(1/6)*(X()+-1im*Y())
basis[13]=sqrt(1/2)*(X()+1im*Y())
basis[14]=-sqrt(1/3)*Z(false)+sqrt(1/3)*(X()+-1im*Y())

#Check orthonormality
#=
orthonormality=zeros(Float64,14,14)
for i=1:14
    for j=1:14
        orthonormality[i,j]=basis[i]*basis[j]
    end
end
println(all(a->abs(a)<1e-2,orthonormality-I))=#

#Define px based on the ket's type
function px14(ket::Ket,P0,P0p,Q)
    if ket.type=='X'
        return 1im*P0*ket.factor*s(ket.spin)
    elseif ket.type=='x'
        return 1im*P0p*ket.factor*s(ket.spin)
    elseif ket.type=='s'
        return -1im*P0*ket.factor*X(ket.spin)+-1im*P0p*ket.factor*x(ket.spin)
    elseif ket.type=='Z'
        return -1im*Q*ket.factor*y(ket.spin)
    elseif ket.type=='y'
        return 1im*Q*ket.factor*Z(ket.spin)
    elseif ket.type=='z'
        return 1im*Q*ket.factor*Y(ket.spin)
    elseif ket.type=='Y'
        return -1im*Q*ket.factor*z(ket.spin)
    end
end
#Make px a linear operator over terms
px14(terms::Terms,args...)=sum(px14.(terms.terms,args...))

function py14(ket::Ket,P0,P0p,Q)
    if ket.type=='X'
        return -1im*Q*ket.factor*z(ket.spin)
    elseif ket.type=='x'
        return 1im*Q*ket.factor*Z(ket.spin)
    elseif ket.type=='s'
        return -1im*P0*ket.factor*Y(ket.spin)+-1im*P0p*ket.factor*y(ket.spin)
    elseif ket.type=='Z'
        return -1im*Q*ket.factor*x(ket.spin)
    elseif ket.type=='y'
        return 1im*P0p*ket.factor*s(ket.spin)
    elseif ket.type=='z'
        return 1im*Q*ket.factor*X(ket.spin)
    elseif ket.type=='Y'
        return 1im*P0*ket.factor*s(ket.spin)
    end
end
py14(terms::Terms,args...)=sum(py14.(terms.terms,args...))

function pz14(ket::Ket,P0,P0p,Q)
    if ket.type=='X'
        return -1im*Q*ket.factor*y(ket.spin)
    elseif ket.type=='x'
        return 1im*Q*ket.factor*Y(ket.spin)
    elseif ket.type=='s'
        return -1im*P0*ket.factor*Z(ket.spin)+-1im*P0p*ket.factor*z(ket.spin)
    elseif ket.type=='Z'
        return 1im*P0*ket.factor*s(ket.spin)
    elseif ket.type=='y'
        return 1im*Q*ket.factor*X(ket.spin)
    elseif ket.type=='z'
        return 1im*P0p*ket.factor*s(ket.spin)
    elseif ket.type=='Y'
        return -1im*Q*ket.factor*x(ket.spin)
    end
end
pz14(terms::Terms,args...)=sum(pz14.(terms.terms,args...))

#Returns pmatx, pmaty, pmatz for given P0, P0p, Q
function pMatrices14(P0,P0p,Q)
    [[basis15[i]*px14(basis15[j],P0,P0p,Q) for j=1:14] for i=1:14], [[basis15[i]*py14(basis15[j],P0,P0p,Q) for j=1:14] for i=1:14], [[basis15[i]*pz14(basis15[j],P0,P0p,Q) for j=1:14] for i=1:14]
end

Xd,Yd,Zd,xc,yc,zc,Xv,Yv,Zv,Sv,Sc,Su,Sq,Dx,Dz=genBasis(["Xd","Yd","Zd","xc","yc","zc","Xᵥ","Yᵥ","Zᵥ","Sᵥ","Sc","Sᵤ","Sq","Dₓ","Dz"])

basis30=zeros(Terms,30)
basis30[1]=Terms(1im*Sq())
basis30[2]=Terms(1im*Sq(false))
basis30[3]=sqrt(1/2)*(Xd()+1im*Yd())
basis30[4]=-sqrt(1/2)*(Xd(false)+-1im*Yd(false))
basis30[5]=-sqrt(2/3)*Zd()+sqrt(1/6)*(Xd(false)+1im*Yd(false))
basis30[6]=sqrt(2/3)*Zd(false)+sqrt(1/6)*(Xd()+-1im*Yd())
basis30[7]=sqrt(1/3)*Zd()+sqrt(1/3)*(Xd()+1im*Yd())*sdown
basis30[8]=-sqrt(1/3)*Zd(false)+sqrt(1/3)*(Xd()+-1im*Yd())
basis30[9]=Terms(1im*Dz())
basis30[10]=Terms(1im*Dz(false))
basis30[11]=Terms(1im*Dx())
basis30[12]=Terms(1im*Dx(false))
basis30[13]=Terms(1im*Su())
basis30[14]=Terms(1im*Su(false))
basis30[15]=sqrt(1/2)*(xc()+1im*yc())
basis30[16]=-sqrt(1/2)*(xc()+-1im*yc())*sdown
basis30[17]=-sqrt(2/3)*zc()+sqrt(1/6)*(xc()+1im*yc())*sdown
basis30[18]=sqrt(2/3)*zc(false)+sqrt(1/6)*(xc()+-1im*yc())
basis30[19]=sqrt(1/3)*zc()+sqrt(1/3)*(xc()+1im*yc())*sdown
basis30[20]=-sqrt(1/3)*zc(false)+sqrt(1/3)*(xc()+-1im*yc())
basis30[21]=Terms(1im*Sc())
basis30[22]=Terms(1im*Sc(false))
basis30[23]=sqrt(1/2)*(Xv()+1im*Yv())
basis30[24]=-sqrt(1/2)*(Xv()+-1im*Yv())*sdown
basis30[25]=-sqrt(2/3)*Zv()+sqrt(1/6)*(Xv()+1im*Yv())*sdown
basis30[26]=sqrt(2/3)*Zv(false)+sqrt(1/6)*(Xv()+-1im*Yv())
basis30[27]=sqrt(1/3)*Zv()+sqrt(1/3)*(Xv()+1im*Yv())*sdown
basis30[28]=-sqrt(1/3)*Zv(false)+sqrt(1/3)*(Xv()+-1im*Yv())
basis30[29]=Terms(1im*Sv())
basis30[30]=Terms(1im*Sv(false))

#Check orthonormality
#=
orthonormality=zeros(ComplexF64,30,30)
for i=1:30
    for j=1:30
        orthonormality[i,j]=basis30[i]*basis30[j]
        if i!=j && abs(orthonormality[i,j])>1e-2
            println((i,j))
        end
    end
end
println(all(a->abs(a)<1e-2,orthonormality-I))=#

function px30(ket::Ket,P0,Pd,P2,P2d,P3,P3d,Ps,Qvc,Pu,Qcd)
    if ket.type=="Xᵥ"
        return 1im*ket.factor*(P0*Sc(ket.spin)+P2*Sq(ket.spin)+P3*Dz(ket.spin))
    elseif ket.type=="Xd"
        return 1im*ket.factor*(Pd*Sc(ket.spin)+P2d*Sq(ket.spin)+P3d*Dz(ket.spin))
    elseif ket.type=="xc"
        return 1im*ket.factor*(Ps*Sv(ket.spin)+Pu*Su(ket.spin))
    else
        return Ket(1,0,true)
    end
end
px30(terms::Terms,args...)=sum(px30.(terms.terms,args...))

function py30(ket::Ket,P0,Pd,P2,P2d,P3,P3d,Ps,Qvc,Pu,Qcd)
    if ket.type=="zc"
        return 1im*ket.factor*Qvc*Xv(ket.spin)
    elseif ket.type=="Zd"
        return 1im*Qcd*ket.factor*xc()
    else
        return Ket(1,0,true)
    end
end
py30(terms::Terms,args...)=sum(py30.(terms.terms,args...))

function pMatrices30(P0,Pd,P2,P2d,P3,P3d,Ps,Qvc,Pu,Qcd)
    [[basis30[i]*px30(basis30[j],P0,Pd,P2,P2d,P3,P3d,Ps,Qvc,Pu,Qcd) for j=1:30] for i=1:30], [[basis30[i]*py30(basis30[j],P0,Pd,P2,P2d,P3,P3d,Ps,Qvc,Pu,Qcd) for j=1:30] for i=1:30]
end

px,py=pMatrices30(9.232,0.195,4.891,9.392,4.328,5.819,3.045,7.998,8.648,4.068)
