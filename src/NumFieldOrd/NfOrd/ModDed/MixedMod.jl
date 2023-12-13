#export MixedMod, quo_pmat, torsion, torsion_free, quo_mixed_mod, ModDedHom, ModDedHom_tfree, image, kernel, inner_sum, direct_sum, product
#TODO rename MixedMod to ModDed and let it include torsion free

###############################################################################
# Submodule of Torsion-free Module
###############################################################################

# checks if a vector is in a torsion free module
function _in_span_(v::AbstractAlgebra.Generic.MatSpaceElem{nf_elem},  P::Hecke.PMat)
   if !(nrows(v)==1) && !(ncols(v)==1)
        return "Please enter a row or column matrix in the first argument"
   end
    @assert length(v) == ncols(P)
   if !can_solve(transpose(P.matrix),transpose(v))
        return false
    else 
        sol = can_solve_with_solution(transpose(P.matrix),transpose(v))[2]
        for i = 1:nrows(P.matrix)
            ideal = P.coeffs[i]
            if !(sol[i,1] in ideal)
                return false
            end
        end
    end
    return true
end

function _in_span_(v::AbstractAlgebra.Generic.MatSpaceElem{nf_elem}, a::Hecke.NfAbsOrdFracIdl{AnticNumberField, nf_elem}, P::Hecke.PMat)
   if !(nrows(v)==1) && !(ncols(v)==1)
        return "Please enter a row or column matrix in the first argument"
   end
    @assert length(v) == ncols(P)
   if !can_solve(transpose(P.matrix),transpose(v))
        return false
    else 
        sol = can_solve_with_solution(transpose(P.matrix),transpose(v))[2]
        for i = 1:nrows(P)
            ideal = P.coeffs[i]*inv(a)
            if !(sol[i,1] in ideal)
                return false
            end
        end
    end
    return true
end



@doc Markdown.doc"""
   spans_subset_of_pseudohnf(M::Hecke.PMat, N::Hecke.PMat)
Checks if M is a submodule of N.
"""
function spans_subset_of_pseudohnf(M::Hecke.PMat, P::Hecke.PMat)
    if nrows(M.matrix)>nrows(P.matrix)
        return false
    end
  for j in 1:nrows(M)
    vv = sub(M.matrix, j:j, 1:ncols(M))
            if !_in_span_(vv,M.coeffs[j],P)
                return false
            end
  end
  return true
end

    
function trans_mat(M::Hecke.PMat, N::Hecke.PMat)
    if !(spans_subset_of_pseudohnf(N,M))
        return "Not a submodule"
    end
    M1 = pseudo_hnf(M)
    N1 = pseudo_hnf(N)
    m1 = nrows(M.matrix)
    m2 = ncols(M.matrix)
    n1 = nrows(N1.matrix)
    n2 = ncols(N1.matrix)
    p = 0
    rr = zeros(Int8,n1,1)
    r = 0
    if !(n2==m2)
        return error("different number of columns")
    end
    s = matrix(base_ring(N1.matrix),n1,m1,[0 for i=1:n1*m1])
    sol = matrix(base_ring(N1.matrix),m2,1,[0 for i=1:m2])
    c = matrix(base_ring(N1.matrix),m2,1,[0 for i=1:m2])
    Mt = transpose(M1.matrix)
    #Mt * sol = c
    for i = 1:n1
        for j in 1:n2
            c[j,1] = N1.matrix[i,j]
        end
        sol = can_solve_with_solution(Mt,c)[2]
        for j in 1:m1
            s[i,j] = transpose(sol)[1,j] 
        end
    end
    for i in 1:n1
        if s[i,1:m1]!=matrix(base_ring(N1.matrix),1,m1,[0 for i in 1:m1])
            p+=1
            rr[i,1] = 1
        end
    end
    s1 = matrix(base_ring(N1.matrix),p,m1,[0 for i=1:p*m1])
    for i in 1:n1
        if rr[i,1]==1
            r+=1
            s1[r,1:m1]=s[i,1:m1]
        end
    end
    return s1
end

#Cohen's algorithm for snf of bimatrix
function cohen_snf_with_transform(P::PMat2) 
    X = deepcopy(P)
    #step 1 - initialise
    unfinished = 1
    step2, step3, step4, step5, step6, step7, step8, step9, step10 = 1,0,0,0,0,0,0,0,0
    a,b = k(), k()
    n = nrows(X.matrix)
    RI = deepcopy(X.row_coeffs)
    CI = deepcopy(X.col_coeffs)
    V = matrix(base_ring(X.matrix), n, n, [0 for i in 1:n*n]) # for row operations
    U = matrix(base_ring(X.matrix), n, n, [0 for i in 1:n*n]) # for column operations
    for i in 1:n
        V[i,i] = 1
        U[i,i] = 1
    end
    rounds = 0
    c = 0
    j = 0
    n == 1 && return (X,U,V)
    i = n 
    while unfinished == 1 && rounds < 6
        while step2 ==1
            j = i
            c = 0
            step2, step3 = 0,1
        end
        while step3 == 1
            if j == 1
                step3, step5 = 0,1
            else
                j = j-1
                if X.matrix[i,j] != 0 
                    step3, step4 = 0,1
                end
            end
        end
        while step4 == 1
            delta = simplify((X.matrix[i,i]*CI[i]) + (X.matrix[i,j]*CI[j]))
            if X.matrix[i,i] == 0 
                u,v = 0, inv(X.matrix[i,j])
            elseif X.matrix[i,j] == 0
                u,v = inv(X.matrix[i,i]), 0
            else
                e,f = idempotents(numerator(simplify(X.matrix[i,i]*(CI[i]*inv(delta)))),numerator(simplify(X.matrix[i,j]*(CI[j]*inv(delta)))))
                u = e // X.matrix[i,i]
                v = f // X.matrix[i,j]
            end
            for count in 1:n
                a, b = deepcopy(X.matrix[i,i]),deepcopy(X.matrix[i,j]) 
                X.matrix[count,j], X.matrix[count,i] = (a*X.matrix[count,j]) - (b*X.matrix[count,i]), u*X.matrix[count,i] + v*X.matrix[count,j]
                U[count,j], U[count,i] = (a*U[count,j]) - (b*U[count,i]), u*U[count,i] + v*U[count,j]
            end
            CI[j], CI[i] = simplify(CI[i]*CI[j]*inv(delta)), simplify(delta)
            step4, step3 = 0,1
        end
        while step5 == 1
            j = i
            if X.matrix[i,i] != 1 
                for count in 1:n
                    U[count,i] = U[count,i]//X.matrix[i,i]
                end
                CI[i] = simplify(CI[i]*X.matrix[i,i])
                X.matrix[i,i] = 1
            end
            step5, step6 = 0,1
        end
        while step6 == 1
            if j == 1 
                step6, step8 = 0,1
            else 
                j = j-1
                if X.matrix[j,i] != 0 
                    step6, step7 = 0,1
                end
            end
        end
        while step7 == 1
            delta = simplify(inv(RI[i]) + X.matrix[j,i]*inv(RI[j]))
            if X.matrix[j,i] == 0
                u,v = 1,0
            else
                u,f = idempotents(numerator(simplify(inv(RI[i]*delta))), numerator(simplify(X.matrix[j,i]*(inv(RI[j]*delta)))))
                v = f // X.matrix[j,i]
            end
            for count in 1:n
                a = X.matrix[j,i]
                X.matrix[j,count], X.matrix[i,count] = X.matrix[j,count] - a*X.matrix[i,count], u*X.matrix[i,count] + v*X.matrix[j,count]
                V[j,count], V[i,count] = V[j,count] - a*V[i,count], u*V[i,count] + v*V[j,count]
            end
            RI[j], RI[i] = simplify(RI[i]*RI[j]*delta), simplify(inv(delta))
            c = c+1
            step7 = 0
            step6 = 1
        end
        while step8 == 1
            if c == 0 
                step8 = 0
                step9 = 1
            else
                step8 = 0
                step2 = 1
            end
        end
        while step9 == 1
            b = simplify(CI[i]*inv(RI[i]))
            for k in 1:i
                for l in 1:i
                    for gen in gens(simplify(X.matrix[k,l]*CI[l]*inv(RI[k])))
                        if !(gen in b) 
                            delta = RI[i]*inv(RI[k])
                            d = gens(simplify(delta))[1]
                            for count in 1:n
                                X.matrix[i,count] = X.matrix[i,count] + d*X.matrix[k,count]
                                V[i,count] = V[i,count] + d*V[k,count]
                            end      
                            step9, step2 = 0,1
                            rounds+=1
                            @goto escape_step9
                        end
                    end
                end
            end
            step9,step10 = 0,1
        end
        @label escape_step9
        while step10 == 1
            if i >= 3
                i = i-1
                step10,step2 = 0,1
            else
                for count in 1:n
                    U[count,1] = U[count,1]//X.matrix[1,1]
                    CI[1] = CI[1]*X.matrix[1,1]
                    X.matrix[1,1] = 1
                    unfinished = 0
                    return V,PMat2(X.matrix, RI, CI),U
                end
            end
        end
    end
    return "unfinished"
end

function cohen_snf(P::PMat2)
    for i in 1:length(cohen_snf_with_transform(P))
        if typeof(cohen_snf_with_transform(P)[i])==PMat2
            return cohen_snf_with_transform(P)[i]
        end
    end
end
    

function quot_pmat2_snf(M::Hecke.PMat, N::Hecke.PMat)
    M1 = pseudo_hnf(M)
    N1 = pseudo_hnf(N)
    n = nrows(N.matrix)
    m = nrows(M.matrix)
    a = deepcopy(pseudo_hnf(M).coeffs)
    b = deepcopy(pseudo_hnf(N).coeffs)
    T = trans_mat(M,N)
    @assert pseudo_hnf(pseudo_matrix(T, N1.coeffs)) == pseudo_matrix(T, N1.coeffs) # show that this should be true AND that this suffices
    rnk = rank(T)
    TT = PMat2(T[1:rnk, 1:rnk], N1.coeffs[1:rnk], M1.coeffs[1:rnk])
    return cohen_snf(TT)
end

###############################################################################
# Mixed Module
###############################################################################

#Quotient module M/N for torsion free modules M and N
function mixed_mod_snf(M::Hecke.PMat, N::Hecke.PMat)
    M1 = pseudo_hnf(M)
    N1 = pseudo_hnf(N)
    rnk = rank(trans_mat(M,N))
    T1 = pseudo_matrix(matrix(base_ring(M.matrix), nrows(M) - rnk, ncols(M), [0 for i in 1:(nrows(M) - rnk)*ncols(M)]), M1. coeffs[rnk+1:nrows(M)])
    for i in 1:nrows(M) - rnk
        T1.matrix[i,i] = 1 
    end
    T2 = Vector{NfOrdQuoRing}(undef, 0)
    X = quot_pmat2_snf(M,N) 
    for i in 1:nrows(X.matrix)
        ideali = simplify(X.row_coeffs[i]*inv(X.col_coeffs[i])) 
        @assert denominator(ideali) == 1
        ideali != ideal(base_ring(M),1) && push!(T2, ResidueRing(base_ring(M), numerator(ideali)))
    end
    return T1,T2
end    
    


@doc Markdown.doc"""
    MixedMod
The type of mixed modules over Dedekind domains.
Mixed modules are determined as a quotient of torsion-free modules (pseudomatrices)
and by their base ring. 
Mixed Modules have a torsion free submodule and a purely torsion submodule.
"""

mutable struct MixedMod #maybe rename to ModDed? TODO check for type ModDed
    quo_pmat::Tuple{Hecke.PMat,Hecke.PMat}
    base_ring::NfOrd
    torsion::Vector{NfOrdQuoRing}
    torsion_free::Hecke.PMat
    
    function MixedMod(M::Hecke.PMat, N::Hecke.PMat)
        !(spans_subset_of_pseudohnf(N,M)) && error("not a submodule")
        T = mixed_mod_snf(M,N) #was mixed_mod!
        z = new()
        z.quo_pmat = (M,N)
        z.base_ring = base_ring(M)
        z.torsion = T[2]
        z.torsion_free = T[1]
        return z
    end
end

function mixed_mod(M::Hecke.PMat, N::Hecke.PMat)
    return MixedMod(M,N)
end

function base_ring(P::MixedMod)
    return P.base_ring
end
function quo_pmat(P::MixedMod)
    return P.quo_pmat
end
function torsion(P::MixedMod)
    return P.torsion
end
function torsion_free(P::MixedMod)
    return P.torsion_free
end

function torsion_dim(P::MixedMod)
    length(P.torsion)
end

function Base.show(io::IO, P::MixedMod)
    if nrows(P.torsion_free.matrix)==0
        print(io, "Torsion Module over the base ring ", "\n", base_ring(P), "\n", "\n")
        if length(P.torsion) == 0
            print(io, "Trivial module")
        else
        print(io, length(P.torsion), "-direct sum of quotient rings")
        end
        for i in 1:length(torsion(P))
            print(io, "\n", "Base ring quotiented by the ideal <", torsion(P)[i].ideal.gen_one, ", ", torsion(P)[i].ideal.gen_two, ">")
        end
    else
        print(io, "Mixed Module over the base ring ","\n", base_ring(P),"\n")
        print(io, "- - -", "\n", "Torsion free part:", "\n")
        print(io, P.torsion_free,"\n", "- - -", "\n")
        print(io, "and torsion part:", "\n")
        print(io, length(P.torsion), "-direct sum of quotient rings")
        for i in 1:length(P.torsion)
            x = simplify(torsion(P)[i].ideal)
            print(io, "\n", "Base ring quotiented by the ideal <", x.gen_one, ", ", x.gen_two, ">")
        end
    end
end

@doc Markdown.doc"""
    spans_subset_of_mixedmod(N::Hecke.PMat, M::Hecke.PMat)
Checks if N is a submodule of M.
"""
function spans_subset_of_mixedmod(N::MixedMod, M::MixedMod) 
    T1 = quo_pmat(M)
    T2 = quo_pmat(N)
    if T2[2]!=T1[2]
        return false
    elseif !(spans_subset_of_pseudohnf(T2[1],T1[1]))
        return false
    end
    return true
end

@doc Markdown.doc"""
    quo_mixedmod(M::Hecke.PMat, N::Hecke.PMat)
Returns the mixed module which is the quotient module of M with N, where N is a submodule of M.
"""
function quo_mixedmod(M::MixedMod, N::MixedMod)
    !spans_subset_of_mixedmod(N,M) && error("not a submodule")
    T1 = quo_pmat(M)
    T2 = quo_pmat(N)
    return MixedMod(T1[1],T2[1])
end

###################################################
# construction for maps between torsion-free modules
###################################################

#called in constructor for map
#right multiplication of matrix with pmat
function image(P::Hecke.PMat, f::AbstractAlgebra.Generic.MatSpaceElem)
    base_ring(P.matrix) == base_ring(f) || @error "base ring must be same"
    N = pseudo_hnf(pseudo_matrix(P.matrix*f,P.coeffs))
    return N
end


@doc Markdown.doc"""
    ModDedHom_tfree
Homomorphism between torsion free modules over dedekind domain.
Determined by a matrix over the base field.
"""
###this should be a subtype of an abstract type? 
mutable struct ModDedHom_tfree
    domain::Hecke.PMat
    pbasis::Vector{PseudoElem}
    codomain::Hecke.PMat
    matrix::AbstractAlgebra.Generic.MatSpaceElem
    #matrix right multiplication
    
    function ModDedHom_tfree(d::Hecke.PMat, c::Hecke.PMat, f::AbstractAlgebra.Generic.MatSpaceElem, check::Bool = true)
        if check
            N = image(d,f)
            spans_subset_of_pseudohnf(N,c) || @error "codomain error"
        end
        z = new()
        z.domain = d
        z.codomain = c
        z.matrix = f
        z.pbasis = pseudo_basis(d)
        return z
    end
end

function Base.show(io::IO, f::ModDedHom_tfree, showrange::Bool = false)
    print(io, "Module Map over Dedekind Domain", base_ring(f.domain.matrix), "\n")
    print(io, "- - - -", "\n", "with domain:", "\n")
    print(io, "Module generated by pseudo basis", "\n")
    for i in 1:length(f.pbasis)
        print(io, f.pbasis[i], "\n")
    end
    if showrange
        print(io, "- - - -", "\n", "and codomain", "\n")
        print(io, "Module represented by", "\n", f.codomain, "\n")
    end
    print(io, "- - - -", "\n")
    print(io, "defined by the matrix", "\n")
    print(io, f.matrix)
end

@doc Markdown.doc"""
     image(M::ModDedHom_tfree)
Returns the image of f
"""
function image(f::ModDedHom_tfree)
    return image(pseudo_hnf(f.domain), f.matrix)
end


@doc Markdown.doc"""
     kernel(M::ModDedHom_tfree)
Returns the kernel of the map f
"""
# Cohen page 35 point (5). Unsure if the alg already implemented requires reverse order of rows. CHECK. 
function kernel(f::ModDedHom_tfree) 
    R = base_ring(f.domain)
    H, U = pseudo_hnf_with_transform(image(f))
    rows = nrows(H.matrix) - rank(H.matrix)
    if rows == 0 
        return pseudo_matrix(matrix(base_ring(f.domain), 0, 0, []))
    else
        ker = matrix(R, rows , ncols(U),[0 for i in 1:ncols(U)*rows])
        for i in 1:rows
            ker[i, 1:ncols(U)] = U[i, 1:ncols(U)]
        end
        return pseudo_matrix(ker, H.coeffs[1:rows])
    end
end


#to check type of map
function same_module(M::Hecke.PMat, N::Hecke.PMat)
    spans_subset_of_pseudohnf(N,M) && spans_subset_of_pseudohnf(M,N) || return false
    return true
end

function is_surjective(f::ModDedHom_tfree)
    same_module(image(f), f.codomain) || return false
    return true
end

function is_injective(f::ModDedHom_tfree)
    kernel(f) == 0 || return false
    return true
end

function is_isomorphic(f::ModDedHom_tfree)
    same_module(image(f), f.domain) || return false
    return true
end



###################################################
# constructions for maps between mixed modules
###################################################
function pseudo_basis(M::MixedMod)
    n1 = nrows(M.quo_pmat[1])
    n2 = nrows(M.quo_pmat[2])
    H1 = pseudo_hnf(M.quo_pmat[1])
    H2 = pseudo_hnf(M.quo_pmat[2])
    X = Vector{PseudoElem}(undef,n1)
    Y = Vector{PseudoElem}(undef,n2)
    for i in 1:n1
        X[i] = PseudoElem(H1.coeffs[i], [H1.matrix[i,j] for j in 1:ncols(H1)])
    end
    for i in 1:n2
        Y[i] = PseudoElem(H2.coeffs[i], [H2.matrix[i,j] for j in 1:ncols(H2)])
    end
    return [X,Y]
end

@doc Markdown.doc"""
    ModDedHom
Homomorphism between mixed modules over dedekind domain.
Determined by a matrix over the base field.
"""
mutable struct ModDedHom
    domain::MixedMod
    pbasis::Vector{Vector{PseudoElem}}
    codomain::MixedMod
    matrix::AbstractAlgebra.Generic.MatSpaceElem
    #the matrix is the same for pbasis of A and B, where domain = A/B
    
    function ModDedHom(d::MixedMod, c::MixedMod, f::AbstractAlgebra.Generic.MatSpaceElem, check::Bool = true)
        #f: d = A/B -----> c = V/W ; 
            if check
            #the associated ModDedHom_tfree should be injective to avoid 0 division error
            N1 = image(d.quo_pmat[1],f)
            N2 = image(d.quo_pmat[2],f)
            spans_subset_of_pseudohnf(N2, c.quo_pmat[2]) || @error "not well defined" #checks if im(B) < W
            spans_subset_of_pseudohnf(N1, c.quo_pmat[1]) || @error "codomain error" #checks if im(A) < V
        end
        z = new()
        z.domain = d
        z.codomain = c
        z.matrix = f
        z.pbasis = pseudo_basis(d)
        return z 
    end
end

function Base.show(io::IO, f::ModDedHom)
        print(io, "Module Map over Dedekind Domain", base_ring(f.domain.quo_pmat[1].matrix), "\n")
    print(io, "- - - -", "\n", "with domain:", "\n")
    print(io, "Mixed Module", "\n")
    show(io, f.domain)
    print(io, "\n", "- - - -", "\n", "and codomain", "\n")
    print(io, "Mixed Module", "\n")
    show(io, f.codomain)
    print(io, "\n", "- - - -", "\n")
    print(io, "defined by the matrix", "\n")
    print(io, f.matrix)
end

@doc Markdown.doc"""
     pullback(M::ModDedHom_tfree)
Returns the pullback of (the image of) a morphism if it is injective, nothing if not injective.
"""
function pull_back(f::ModDedHom_tfree, M::Hecke.PMat)
    M1 = matrix(base_ring(M), 2*nrows(M), ncols(M),[0 for i in 1:2*nrows(M)*ncols(M)])
    count = 0
    for i in 1:nrows(M)
        for gen in gens(M.coeffs[i])
            count+=1
            M1[count, 1:ncols(M)] = gen*M.matrix[i, 1:ncols(M)]
        end
    end
    V = matrix(base_ring(M),0,ncols(M),[])
    for row in 1:nrows(M1)
        x = can_solve_with_solution(f.matrix, transpose(M.matrix[row,1:ncols(M)]))[2]
        if _in_span_(transpose(x), f.domain)
            V = [V;x]
        end
    end
    return reduce_hnf(pseudo_matrix(V))    
end

function image(f::ModDedHom)
    return MixedMod(image(f.domain.quo_pmat[1], f.matrix), image(f.domain.quo_pmat[2], f.matrix))
end

function kernel(f::ModDedHom) 
    F = ModDedHom_tfree(f.domain.quo_pmat[1], f.codomain.quo_pmat[1], f.matrix)
    return pull_back(F, f.codomain.quo_pmat[2])
end

function same_module(M::MixedMod, N::MixedMod)
    M1, M2 = quot_pmat(M)[1], quot_pmat(M)[2]
    N1, N2 = quot_pmat(N)[1], quot_pmat(N)[2]
    MX, MY = mixed_mod_snf(M1,M2)
    NX, NY = mixed_mod_snf(N1,N2)
    if MY==NY && same_module(MX,NX)
        return true
    end
    return false
end

function is_surjective(f::ModDedHom)
    same_module(image(f), f.codomain) && return true
    return false
end

function is_injective(f::ModDedHom) #must check what it returns
    if kernel(f).torsion == NfOrdQuoRing[] && nrows(kernel(f).torsion_free)==0
        return true
    else
        return false
    end
end

function is_isomorphic(f::ModDedHom)
    if is_injective(f) && is_surjective(f) 
        return true
    end
    return false
end


###################################################
# extra functions
###################################################


function reduce_hnf(m::Hecke.PMat)
    P = deepcopy(pseudo_hnf(m))
    k = 0
    j = 0
    x = zeros(Int8, nrows(m),1)
    for i in 1:nrows(m)
        if P.matrix[i,1:ncols(m)] != matrix(base_ring(m.matrix),1,ncols(m),[0 for i in 1:ncols(m)]) && P.coeffs[i]!=ideal(base_ring(P),0)
            k+=1
            x[i,1] = 1
        end
    end
    W = pseudo_matrix(matrix(base_ring(P.matrix), k, ncols(P), [0 for i in 1:ncols(P)*k]))
    for i in 1:nrows(m)
        if x[i,1]==1
            j+=1
            W.matrix[j,1:ncols(P)] = P.matrix[i,1:ncols(P)]
            W.coeffs[j] = P.coeffs[i]
        end
    end
    return pseudo_hnf(W)
end

#ambient vector space
function ambient_vector_space(M::Hecke.PMat)
    return VectorSpace(base_ring(M.matrix), ncols(M))
end


@doc Markdown.doc"""
     inner_sum(M::Hecke.PMat, N::Hecke.PMat)
Returns module M + N.
"""
function inner_sum(M::Hecke.PMat, N::Hecke.PMat)
    ambient_vector_space(M) == ambient_vector_space(N) || return error("Not submodules of the same ambient space")
    X = pseudo_matrix(matrix(base_ring(M.matrix), nrows(M)+nrows(N), ncols(M), [0 for i in 1:(nrows(M)+nrows(N))*ncols(M)]))
    for i in 1:nrows(N)
        X.coeffs[i] = N.coeffs[i]
        X.matrix[i,1:ncols(N)] = N.matrix[i,1:ncols(N)]
    end
    for i in 1:nrows(M)
        X.coeffs[i+nrows(N)] = M.coeffs[i]
        X.matrix[i+nrows(N), 1:ncols(M)] = M.matrix[i,1:ncols(M)]
    end
    return reduce_hnf(X)
end

@doc Markdown.doc"""
     direct_sum(M::Hecke.PMat, N::Hecke.PMat)
Returns a pseudo matrix representing the direct sum of the modules represented by M and N.
"""
function direct_sum(M::Hecke.PMat, N::Hecke.PMat)
    base_ring(M) == base_ring(N) || return error("base ring must be same")
    X = matrix(base_ring(M.matrix), nrows(M)+nrows(N), ncols(M)+ncols(N), [0 for i in 1:(nrows(M)+nrows(N))*(ncols(M)+ncols(N))])
    I = Vector{Hecke.NfAbsOrdFracIdl}(undef, nrows(M)+nrows(N))
    for i in 1:nrows(M)
        for j in 1:ncols(M)
            X[i,j] = M.matrix[i,j]
            I[i] = M.coeffs[i]
        end
    end
    for i in 1:nrows(N)
        for j in 1:ncols(N)
            X[nrows(M)+i, ncols(N)+j] = N.matrix[i,j]
            I[nrows(M)+i] = N.coeffs[i]
        end
    end
    return pseudo_matrix(X,I)
end



@doc Markdown.doc"""
     product(M::Hecke.PMat, N::Hecke.PMat)
Returns product of modules represented by M and N.
"""
function product(M::Hecke.PMat, N::Hecke.PMat)
    if ambient_vector_space(M) != ambient_vector_space(N)
        return error("incompatible")
    end
    M1 = pseudo_hnf(M)
    N1 = pseudo_hnf(N)
    x = 0
    V = Vector{Hecke.NfOrdFracIdl}(undef, nrows(M)*nrows(N))
    X = matrix(base_ring(M.matrix), nrows(M)*nrows(N), ncols(M), [0 for i in 1:nrows(M)*nrows(N)*ncols(M)])
    for i in 1:nrows(M)
        for j in 1:nrows(N)
            x+=1
            for t in 1:ncols(M)
                X[x,t] = M.matrix[i,t]*N.matrix[j,t]
            end
            V[x] = M.coeffs[i]*N.coeffs[j]
        end
    end
    return reduce_hnf(pseudo_matrix(X,V))
end


#@doc Markdown.doc"""
#     intersection(M::Hecke.PMat, N::Hecke.PMat)
#Returns intersection of modules represented by M and N.
#"""

#Cohen 1.5.2 point (6)
#broken
function intersection(M::Hecke.PMat, N::Hecke.PMat) #alg from Cohen pg36
    (base_ring(M) == base_ring(N) && nrows(M)==nrows(N) && ncols(M)==ncols(N)) || return error("base ring, rows and columns must be same")
    n = nrows(M)
    m = ncols(M)
    X = matrix(base_ring(M.matrix), 2*n, 2*m, [0 for i in 1:4*n*m])
    I = Vector{Hecke.NfAbsOrdFracIdl{AnticNumberField, nf_elem}}(undef, 2*n)
    X[1:n, 1:m] = M.matrix
    X[n+1:2*n,1:m] = M.matrix
    X[n+1:2*n, m+1:2*m] = N.matrix
    I[1:n] = M.coeffs
    I[n+1:2*n] = N.coeffs
    H = pseudo_hnf(pseudo_matrix(X,I))
    return pseudo_matrix(H.matrix[1:n, 1:m], H.coeffs[1:n])
end