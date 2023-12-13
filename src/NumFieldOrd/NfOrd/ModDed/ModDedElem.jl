export PseudoElem, MixedModElem

###############################################################################
# Elements of Modules over Dedekind Domains
###############################################################################


###############################################################################
# Element of torsion-free module 
###############################################################################

@doc Markdown.doc"""
    PseudoElem
The type of element of a torsion-free module over a Dedekind domain.
Pseudo elements are determined by their base ring, their ideal, and a vector of `nf_elem`. 
"""

mutable struct PseudoElem
    base_ring::NfOrd
    ideal::Hecke.NfAbsOrdFracIdl
    vector::Vector{nf_elem}
    dimension::Int64
    
    function pseudo_elem(I::Union{Hecke.NfAbsOrdFracIdl{AnticNumberField, nf_elem},NfOrdIdl}, v::Vector{nf_elem})
        length(v)==0 && error
        if typeof(I)==NfOrdIdl
            I = fractional_ideal(I)
        end
        if parent(v[1])!==nf(order(I))
            return error
        end
        n = size(v)[1]
        for i in 1:n
            if parent(v[1])!==parent(v[i])
                return error
            end
        end
        z = new()
        z.ideal = I
        z.vector = v
        z.dimension = length(v)
        z.base_ring = order(I)
        return z
    end
    
    function PseudoElem(I::Union{Hecke.NfAbsOrdFracIdl{AnticNumberField, nf_elem},NfOrdIdl}, v::Vector{nf_elem})
        length(v)==0 && error
        if typeof(I)==NfOrdIdl
            I = fractional_ideal(I)
        end
        if parent(v[1])!==nf(order(I))
            return error
        end
        n = size(v)[1]
        for i in 1:n
            if parent(v[1])!==parent(v[i])
                return error
            end
        end
        z = new()
        z.ideal = I
        z.vector = v
        z.dimension = length(v)
        z.base_ring = order(I)
        return z
    end
end
    
function dimension(x::PseudoElem)
    return x.dimension
end
function ideal(x::PseudoElem)
    return x.ideal
end
function vector(x::PseudoElem)
    return x.vector
end
function base_ring(x::PseudoElem)
    return x.base_ring
end

function Base.show(io::IO, e::PseudoElem)
    print(io, "Pseudo Element over ", e.base_ring, "\n")
    print(io, "with ideal 1//", denominator(e.ideal), " * ")
    print(io, numerator(e.ideal), "\n")
    print(io, "and vector [")
    for i in 1:e.dimension
        if i != e.dimension 
            print(io, e.vector[i], ", ")
        else 
            print(io, e.vector[i])
        end
    end
    print(io, "]")
end

@doc Markdown.doc"""
    contained_in(x::PseudoElem, T::Hecke.PMat)
Checks if `x` is contained in `T`.
"""
function contained_in(x::PseudoElem, T::Hecke.PMat)
    if !(base_ring(T)==x.base_ring) 
        return false
    end
    if typeof(x.ideal)==NfOrdIdl
        y=fractional_ideal(x.ideal)
    else
        y=deepcopy(x.ideal)
    end
    X=matrix(base_ring(T.matrix),1,x.dimension,[0 for i in 1:x.dimension])
    for i=1:x.dimension
        X[1,i]=x.vector[i]
    end
    M=pseudo_matrix(X,[y])
    i=_spans_subset_of_pseudohnf_(M,T)
    return i
end

###############################################################################
# Element of Mixed Modules
###############################################################################

@doc Markdown.doc"""
    MixedModElem
The type of element of a mixed module over a Dedekind domain.
Mixed module elements are determined by their order, their cyclic part, which is a vector of `NfOrdQuoRingElem` 
and their non cyclic part, which is a `PseudoElem`. 
"""
mutable struct MixedModElem
    order::NfOrd
    cyclic::Vector{NfOrdQuoRingElem}
    non_cyclic::PseudoElem
    
    function mixed_mod_elem(x::PseudoElem, y::Vector{NfOrdQuoRingElem})
        if length(y)==0
            return x
        end
        for i in 1:length(y)
            if parent(y[1].elem) !=parent(y[i].elem)
                return error("base ring must be same")
            end
        end
        if !(parent(y[1].elem)==base_ring(x))
            return error("base ring must be same")
        end
        z = new()
        z.cyclic = y
        z.non_cyclic = x
        z.order = x.base_ring
        return z
    end
    
    function MixedModElem(x::PseudoElem, y::Vector{NfOrdQuoRingElem})
        if length(y)==0
            return x
        end
        for i in 1:length(y)
            if parent(y[1].elem) !=parent(y[i].elem)
                return error("base ring must be same")
            end
        end
        if !(parent(y[1].elem)==base_ring(x))
            return error("base ring must be same")
        end
        z = new()
        z.cyclic = y
        z.non_cyclic = x
        z.order = x.base_ring
        return z
    end
end

function order(x::MixedModElem)
    return x.order
end

function cyclic(x::MixedModElem)
    return x.cyclic
end

function non_cyclic(x::MixedModElem)
    return x.non_cyclic
end

function Base.show(io::IO, e::MixedModElem)
    print(io, "Mixed Module Element over ", e.order, "\n")
    print(io, "- - -", "\n", "with pseudo element:", "\n")
    print(io, e.non_cyclic)
    print(io, "\n", "- - -", "\n")
    print(io, "and vector of Quotient Ring elements", "\n")
    for i in 1:length(e.cyclic)
        print(io, "(", e.cyclic[i], ") in Residue Ring: ", e.order, " quotiented by ")
        print(io, "<", gens(ideal(parent(e.cyclic[i])))[1], ", ",  gens(ideal(parent(e.cyclic[i])))[2], ">", "\n")
    end
end

@doc Markdown.doc"""
    contained_in(x::MixedModElem, T::MixedMod)
Checks if `x` is contained in `T`.
"""
#the order of the elements has to be the order of the residue fields
function contained_in(x::MixedModElem, T::MixedMod)
    R = base_ring(T.torsion_free)
    if !(contained_in(x.non_cyclic,T.torsion_free)) || length(x.cyclic)!=length(T.torsion)
        return false
    end
    for i=1:length(x.cyclic)
        if parent(x.cyclic[i])!=T.torsion[i]
            return false
        end
    end
    return true
end