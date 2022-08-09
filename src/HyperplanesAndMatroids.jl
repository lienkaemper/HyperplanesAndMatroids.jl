module HyperplanesAndMatroids


#requirements
using Combinatorics #we need this for powerset() function
using LinearAlgebra

__precompile__()

export SignVector, from_hyperplanes, chirotope, sep


#---------------Sign Vector Stuff------------------------#

"""
    SignVector

A SignVector represents a vector w/ real entries,
designed to be equal if and only if they have the same signs

Example:

```julia> SignVector([0, 1, -1, -1])
0+--```



"""
struct SignVector
    vec::Array{Int8}
end

"""
    SignVector(str::String)

Construct a sign vector from a string with characters '0', '+', '-'

Examples:

```julia> SignVector("0+--")
0+--```

"""
function SignVector(str::String)
    vec = zeros(length(str))
    for i in 1:length(str)
        if str[i] == '+'
            vec[i] = 1
        elseif str[i] == '-'
            vec[i] = -1
        end
    end
    return SignVector(vec)
end

"""
    Base.show(io::IO, X::SignVector)

Pretty printing for sign vectors

"""
function Base.show(io::IO, X::SignVector)
    str = ""
    V = X.vec
    for i in V
        if i < 0
            str = str*"-"
        elseif i > 0
            str = str*"+"
        elseif i == 0
            str = str*"0"
        end
    end
    print(io, str)
end

Base.getindex(X::SignVector, i::Int)= X.vec[i]
Base.getindex(X::SignVector, I::Array{Int})= SignVector(X.vec[I])
Base.getindex(X::SignVector, u::UnitRange{Int})= SignVector(X.vec[u])
Base.length(X::SignVector)= length(X.vec)
Base.:-(X::SignVector)= SignVector(-X.vec)
Base.copy(X::SignVector) = SignVector(copy(X.vec))
Base.:(==)(X::SignVector, Y::SignVector) = Int.(sign.(X.vec)) == Int.(sign.(Y.vec))
Base.isequal(X::SignVector, Y::SignVector) = Int.(sign.(X.vec)) == Int.(sign.(Y.vec))
Base.hash(X::SignVector) = hash(Int.(sign.(X.vec)))
Base.last(X::SignVector) = last(X.vec)
Base.:*(rr::Real, X::SignVector) = SignVector(Int.(sign.(rr *X.vec)))

function Base.setindex!(X::SignVector, val::Real, i::Int)
    X.vec[i] = Int(sign(val))
end

function Base.setindex!(X::SignVector, val::Array{Real}, i::Array{Int})
    X.vec[i] = Int(sign(val))
end

function Base.setindex!(X::SignVector, val::SignVector, i::Array{Int})
    X.vec[i] = val.vec
end


#----------------Sign Vector Tools--------------------------#

"""
    positivePart(X::SignVector)

Returns the positive part of a sign vector, i.e. {i ∣ X_i = +}

"""
function positivePart(X::SignVector)
    return(findall(x-> x> 0, X.vec))
end


"""
    negativePart(X::SignVector)

Returns the positive part of a sign vector, i.e. {i ∣ X_i = -}

"""
function negativePart(X::SignVector)
    return(findall(x-> x< 0, X.vec))
end

"""
    support(X::SignVector)

Returns the support of a sign vector, i.e. {i ∣ X_i ≠ 0}

"""
function support(X::SignVector)
    return(findall(x-> x!= 0, X.vec))
end

"""
    zeroPart(X::SignVector)

Returns the zero part of a sign vector, i.e. {i ∣ X_i = 0}

"""
function zeroPart(X::SignVector)
    return(findall(x-> x== 0, X.vec))
end

function orthogonal(X::SignVector, Y ::SignVector)
    x = X.vec
    y = Y.vec
    n = length(x)

    #first, check if supports are disjoint
    if sum([abs(x[i]*y[i]) for i = 1:n]) == 0
        return true
    end

    #otherwise, check that vectors are neither equal nor opposite
    diff = 0
    same = 0
    for i = 1:n
        if x[i]*y[i] == 1
            same = 1
        elseif x[i]*y[i] == -1
            diff = 1
        end
    end
    orth = diff*same
    return Bool(orth)
end




function composition(X::SignVector, Y::SignVector)
    result = copy(X)
    for i = 1:length(X)
        if X[i]==0
            result[i] = Y[i]
        end
    end
    return result
end

function composition(X::SignVector...)
    result = X[1]
    n_inputs = length(X)
    for i = 2:n_inputs
        result = composition(result, X[i])
    end
    return result
end

function sep(X::SignVector,Y::SignVector)
    result = []
    for i =1:length(X)
        if X[i]*Y[i] == -1
            append!(result,i)
        end
    end
    return result
end

#labels a sign vector of length n w/ a number between 0, 2^n-1
function binaryLabel(X::SignVector)
    n = length(X)
    return sum([2^(n-i)*Int64((X[i]+1)/2) for i=1:length(X)])
end

#inverse of binaryLabel
#k: binary label of sign vec
#n: desired length of sign vec
function labelToSignVec(k,n)
    v = last(bitstring(k), n)
    n = length(v)
    result = zeros(Int64, n)
    for i = 1:n
        if v[i] == '0'
            result[i] = -1
        elseif v[i]=='1'
            result[i] = 1
        end
    end
    return SignVector(result)
end



# returns a vector containing all 2^n vectors of length n w/ entries ±1
function allPM(n)
    vecs = [labelToSignVec(i,n) for i =0:2^n-1]
    return vecs
end


#---------------Oriented Matroid Basics----------------------#

mutable struct OrientedMatroid
    topes::Array{SignVector}
    circuits::Dict{Array{Int}, Array{SignVector}}
    cocircuits::Dict{Array{Int}, Array{SignVector}}
    chirotope::Dict{Array{Int}, Int}
    n::Int
    r::Int
    OrientedMatroid() = new()
end

function setGround!(M::OrientedMatroid, n::Int)
    M.n::n
end

function setRank!(M::OrientedMatroid, r::Int)
    M.r::r
end

function setTopes!(M,topes::Vector)
    M.topes = SignVector.(topes)
    return M
end

function setTopes!(M,topes::String)
    M.topes = SignVector.(String.(strip.(split(topes, ','))))
    return M
end

function setCircuits!(M,circuits::Vector)
    values = circuits
    M.circuits = Dict{Array{Int}, Array{SignVector}}()
    for circuit in values
        supp = support(circuit)
        if haskey(M.circuits, supp)
            append!(M.circuits[supp], [circuit])
        else
            M.circuits[supp] = [circuit]
        end
    end
    return M
end

function setCircuits!(M,circuits::String)
    values = SignVector.(String.(strip.(split(circuits, ','))))
    M.circuits = Dict{Array{Int}, Array{SignVector}}()
    for circuit in values
        supp = support(circuit)
        if haskey(M.circuits, supp)
            append!(M.circuits[supp], [circuit])
        else
            M.circuits[supp] = [circuit]
        end
    end
    return M
end

function setCoCircuits!(M,cocircuits::Vector)
    values = cocircuits
    M.cocircuits = Dict{Array{Int}, Array{SignVector}}()
    for cocircuit in values
        supp = support(cocircuit)
        if haskey(M.cocircuits, supp)
            append!(M.cocircuits[supp], [cocircuit])
        else
            M.cocircuits[supp] = [cocircuit]
        end
    end
    return M
end

function setCoCircuits!(M,cocircuits::String)
    values = SignVector.(String.(strip.(split(cocircuits, ','))))
    M.cocircuits = Dict{Array{Int}, Array{SignVector}}()
    for cocircuit in values
        supp = support(cocircuit)
        if haskey(M.cocircuits, supp)
            append!(M.cocircuits[supp], [cocircuit])
        else
            M.cocircuits[supp] = [cocircuit]
        end
    end
    return M
end


#---------------Affine oriented matroids----------------------#

function affineTopes(M; k = M.n)
    aff_topes = []
    actual_planes = setdiff(collect(1:M.n), k)
    for tope in M.topes
        if tope[k] > 0
            push!(aff_topes, tope[actual_planes])
        end
    end
    return aff_topes
end


#---------------Converting Between Representations------------------#

function comp_sign(σ, e)
    r = length(σ)
    current = ((r + 1) %2)*2 -1
    original = (minimum([findall(σ .> e); r+1]) % 2)*2 -1
    return current * original
end

#returns the unique signed cocircuit which is disjoint from χ∖e, positve on e
function basic_cocircuit(χ, σ, e,n)
    signature =  SignVector(zeros(n))
    δ = setdiff(σ , [e])
    a = [δ ;[e]]
    signature[e] = χ[σ] * comp_sign(δ,e)
    for f in setdiff(collect(1:n) ,σ)
        b = [δ ;[f]]
        signature[f] = χ[sort(b)] * comp_sign(δ, f)
    end
    return signature
end

#returns the unique signed cocircuit which is disjoint from χ∖e, positve on e
function basic_circuit(χ, σ, e,n)
    signature =  SignVector(zeros(n))
    signature[e] = 1
    for f in σ
        δ = setdiff(σ, [f])
        a = [δ ;[e]]
        b = [δ ;[f]]
        signature[f] = -1 * χ[sort(b)] * comp_sign(δ, f)* χ[sort(a)] * comp_sign(δ, e)
    end
    return signature
end

#given chirotope, return cocircuits as list of sign vectors
function cocircuits(χ, n)
    cocircs = Vector{SignVector}()
    for (σ, sgn) in χ
        if sgn != 0
            for e in σ
                signature = basic_cocircuit(χ, σ, e, n)
                if signature ∉ cocircs
                    append!(cocircs, [signature, -signature])
                end
            end
        end
    end
    return cocircs
end

#given chirotope, return cocircuits as list of sign vectors
function circuits(χ, n)
    circs = Vector{SignVector}()
    for (σ, sgn) in χ
        if sgn != 0
            for e in setdiff(collect(1:n), σ)
                signature = basic_circuit(χ, σ, e, n)
                if signature ∉ circs
                    append!(circs, [signature, -signature])
                end
            end
        end
    end
    return circs
end

#return topes of rank r matroid on ground set n with chirotope χ.
function topes(χ, n, r)
    topes = []
    for (σ, sgn) ∈ χ
        if sgn != 0
            for α ∈ allPM(r)
                T_σα = composition([α[i] * basic_cocircuit(χ, σ, σ[i], n) for i = 1:r]...)
                if T_σα ∉ topes
                    push!(topes, T_σα)
                end
            end
        end
    end
    return topes
end

# M an oriented matroid with circuits filled in
# fills in topes of M, also returns them
function circuitsToTopes!(M::OrientedMatroid)
    # only need to take one circuit from each opposite pair,
    # since orthogonality X⟂Y iff and only if X⟂-Y.
    circuits = [circ[1] for circ in values(M.circuits)]
    n = length(circuits[1])
    topes = SignVector[]
    # something is a tope if it's orthogonal to all circuits
    for T in allPM(n)
        orth = true
        for C in circuits
            if !orthogonal(C, T)
                orth = false
                break
            end
        end
        if orth == true
            append!(topes, [T])
        end
    end
    M.topes = topes
    return topes
end

# find a list of circuits on support orthogonal to topes
# if topes actually satisfies tope axioms, this will be one pair of opposites
# otherwise, will have more pairs
function findCircuit(topes, support)
    n = length(support)
    m = length(topes[1])
    restriction = [T[support] for T in topes]

    patterns_found = zeros(2^(n))
    count = 0
    i = 1
    while count < 2^(n) && i <= length(restriction)
        label = binaryLabel(restriction[i])
        if patterns_found[label+1]==0
            count += 1
            patterns_found[label+1] = 1
        end
        i+=1
    end
    circuits = SignVector[]
    if count < 2^n
        circuit_labels = findall(iszero, patterns_found)
        for label in circuit_labels
            circuit = SignVector(zeros(m))
            circuit[support] = labelToSignVec(label-1,n)
            append!(circuits,[circuit])
        end
    end
    return circuits
end


# M an oriented matroid with topes filled in
# need to specify a rank, currently only works w/ uniform matroids
# computes circuits from the topes and fills them in
function topesToCircuits!(M; uniform = true)
    rank = M.r
    topes = M.topes
    n = length(topes[1])
    circuits = Dict{Array{Int}, Array{SignVector}}()
    if uniform
        for sigma in powerset(collect(1:n), rank+1, rank+1)
            circuits[sigma] = findCircuit(topes, sigma)
        end
    end
    M.circuits = circuits
    return circuits
end

function chirotopeToTopes!(M)
    χ = M.chirotope
    n = M.n
    r = M.r
    M.topes = topes(χ, n, r)
end

function chirotopeToCoCircuits!(M)
    χ = M.chirotope
    n = M.n
    r = M.r
    setCoCircuits!(M, cocircuits(χ, n))
end

function chirotopeToCircuits!(M)
    χ = M.chirotope
    n = M.n
    r = M.r
    setCircuits!(M, circuits(χ, n))
end

#input: hyperplane arrangement, rows are hyperplane normals
#output: dictionary, σ -> sign(det(ha_sigma))

#this will work, even if matroid is not uniform
function chirotope(ha)
    (n, r) = size(ha)
    chirotope =Dict()
    for σ in powerset(collect(1:n), r, r) #all subsets of 1:n of size r
        σ_2 = sort(σ)
        chirotope[σ_2] = sign(det(ha[σ_2, :]))
    end
    return chirotope
end

function from_hyperplanes(ha)
    M = OrientedMatroid()
    M.n = size(ha, 1)
    M.r = rank(ha)
    M.chirotope = chirotope(ha)
    chirotopeToTopes!(M)
    chirotopeToCoCircuits!(M)
    chirotopeToCircuits!(M)
    return M
end


end
