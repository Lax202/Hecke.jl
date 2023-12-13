# Examples

In this chapter we will construct examples to understand the functionality of modules over dedekind domains. 


```julia
using Oscar
```

     -----    -----    -----      -      -----   
    |     |  |     |  |     |    | |    |     |  
    |     |  |        |         |   |   |     |  
    |     |   -----   |        |     |  |-----   
    |     |        |  |        |-----|  |   |    
    |     |  |     |  |     |  |     |  |    |   
     -----    -----    -----   -     -  -     -  
    
    ...combining (and extending) ANTIC, GAP, Polymake and Singular
    Version[32m 0.11.0 [39m... 
     ... which comes with absolutely no warranty whatsoever
    Type: '?Oscar' for more information
    (c) 2019-2022 by The Oscar Development Team


## Torsion-free modules

First We will construct a Dedekind domain (a maximal order of a number field).


```julia
k,a= quadratic_field(10)
zk = maximal_order(k)
```




    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]



We can then construct torsion-free modules over $zk$. 


```julia
m1=pseudo_matrix(matrix(zk,4,4,[1 0 0 0; 0 2 0 0; 0 0 3 0;0 0 0 1]))
m2=pseudo_matrix(matrix(zk,2,4,[5 0 0 0; 0 6 0 0]))
```




    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <1, 1> with row [5 0 0 0]
    1//1 * <1, 1> with row [0 6 0 0]



We can use the function `spans_subset_of_pseudohnf` to check if m2 is a submodule of m1.


```julia
spans_subset_of_pseudohnf(m2,m1)
```




    true



## Mixed Modules

Then we can take the quotient module of m1 by m2 with the constructor `MixedMod`.


```julia
X = MixedMod(m1,m2)
```




    Mixed Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    - - -
    Torsion free part:
    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <3, 3> with row [1 0 0 0]
    1//1 * <1, 1> with row [0 1 0 0]
    - - -
    and torsion part:
    1-direct sum of quotient rings
    Base ring quotiented by the ideal <15, 165>



This returns a module $X$ of type `MixedMod`, a mixed module over $zk$. Therefore $X$ has a torsion-free submodule (a pseudomatrix):


```julia
X.torsion_free
```




    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <1, 1> with row [1 0 0 0]
    1//1 * <1, 1> with row [0 1 0 0]



and a torsion submodule:


```julia
X.torsion
```




    2-element Vector{NfOrdQuoRing}:
     Quotient of Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
     Quotient of Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]



which is a vector of quotient rings of $zk$.

Crucially, every object of type `MixedMod` has a field which indicates how it was constructed, as the quotient of two torsion-free modules. This field is called `quo_pmat`.


```julia
quo_pmat(X)
```




    (Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <1, 1> with row [1 0 0 0]
    1//1 * <1, 1> with row [0 2 0 0]
    1//1 * <1, 1> with row [0 0 3 0]
    1//1 * <1, 1> with row [0 0 0 1], Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <1, 1> with row [5 0 0 0]
    1//1 * <1, 1> with row [0 6 0 0])





We can define another mixed module that contains $X$.


```julia
m3 = pseudo_matrix(matrix(zk,3,4,[1 0 0 0; 0 2 0 0; 0 0 6 0]))
Y = MixedMod(m3,m2)
```




    Mixed Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    - - -
    Torsion free part:
    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <36, 6> with row [1 0 0 0]
    - - -
    and torsion part:
    1-direct sum of quotient rings
    Base ring quotiented by the ideal <15, 165>




```julia
spans_subset_of_mixedmod(Y,X)
```




    true



Since $Y < X$ as modules, we can create a new module Z which is the quotient of X and Y with the function `quo_mixedmod`


```julia
quo_mixedmod(X,Y)
```




    Mixed Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    - - -
    Torsion free part:
    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <1, 1> with row [1 0 0 0]
    - - -
    and torsion part:
    1-direct sum of quotient rings
    Base ring quotiented by the ideal <2, 2>



A more general example of MixedMod:


```julia
g = pseudo_matrix(matrix(k,4,4,[rand(k, -10:10) for i in 1:16]))
```




    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <1, 1> with row [-10*sqrt(10)-1 -2*sqrt(10)-1 -8*sqrt(10)+4 9*sqrt(10)-7]
    1//1 * <1, 1> with row [-7*sqrt(10) 3*sqrt(10)-8 sqrt(10)+6 7*sqrt(10)+1]
    1//1 * <1, 1> with row [4*sqrt(10)-4 2*sqrt(10)-2 -sqrt(10)-6 sqrt(10)-5]
    1//1 * <1, 1> with row [-sqrt(10)+3 9*sqrt(10)+4 -2*sqrt(10)+10 9*sqrt(10)+8]




```julia
h = pseudo_matrix(matrix(k,4,4,[rand(k, -1:10) for i in 1:16]))
```




    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <1, 1> with row [7*sqrt(10)+6 3*sqrt(10)-1 7*sqrt(10) 2*sqrt(10)+4]
    1//1 * <1, 1> with row [9*sqrt(10)+9 6*sqrt(10) 6*sqrt(10)+4 9*sqrt(10)+9]
    1//1 * <1, 1> with row [sqrt(10)+10 4*sqrt(10)+7 -sqrt(10)+7 2*sqrt(10)+7]
    1//1 * <1, 1> with row [3*sqrt(10)+2 10*sqrt(10)+1 9*sqrt(10)+4 3*sqrt(10)+6]




```julia
MixedMod(inner_sum(g,h), h)
```




    Torsion Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    
    1-direct sum of quotient rings
    Base ring quotiented by the ideal <13908623994, 383521114478567483*sqrt(10) + 167512674736388582758>



There is still an outstanding issue with intersection


```julia
MixedMod(g, intersection(g,h))
```




    Torsion Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    
    Trivial module





## Elements 

An element of a torsion-free module is of the data type `PseudoElement`. A pseudo element is constructed with an ideal of type `NfAbsOrdFracIdl` or `NfOrdIdl` and a vector of elements of the number field `nf_elem`


```julia
x = PseudoElem(ideal(zk,(zk)(2)), [k(), k(2), k(), k(0)])
```




    Pseudo Element over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    with ideal 1//1 * <4, 2>
    Norm: 4
    principal generator 2
    two normal wrt: 4
    and vector [0, 2, 0, 0]




```julia
contained_in(x, m1), contained_in(x,m2)
```




    (true, false)



An element of a mixed module is of the data type `MixedModElem` and is constructed with a pseudo element and a vector of elements of quotient rings `NfQuoRingElem`.


```julia
t = ModDedElem(x, [(ResidueRing(zk,ideal(zk,(zk)(5))))(2),(ResidueRing(zk,ideal(zk,(zk)(3))))(2)])
```




    Mixed Module Element over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    - - -
    with pseudo element:
    Pseudo Element over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    with ideal 1//1 * <4, 2>
    Norm: 4
    principal generator 2
    two normal wrt: 4
    and vector [0, 2, 0, 0]
    - - -
    and vector of Quotient Ring elements
    (2) in Residue Ring: Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)] quotiented by <25, 5>
    (2) in Residue Ring: Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)] quotiented by <9, 3>





```julia
contained_in(t, X), contained_in(t,Y)
```




    (true, false)





## Maps

A map between torsion-free modules is of the data type `ModDedHom_tfree`. It is defined (from a domain to a codomain) by a matrix over the corresponding number field. 


```julia
F = ModDedHom_tfree(m2,m1,matrix(k,4,4,[0 0 3 0; 0 3 0 0; 5 0 0 0; 1 1 1 1]))
```




    Module Map over Dedekind DomainReal quadratic field defined by x^2 - 10
    - - - -
    with domain:
    Module generated by pseudo basis
    PseudoElem(Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)], 1//1 * <25, 5>
    Norm: 25
    principal generator 5
    two normal wrt: 25, nf_elem[1, 0, 0, 0], 4)
    PseudoElem(Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)], 1//1 * <36, 6>
    Norm: 36
    principal generator 6
    two normal wrt: 36, nf_elem[0, 1, 0, 0], 4)
    - - - -
    defined by the matrix
    [0 0 3 0; 0 3 0 0; 5 0 0 0; 1 1 1 1]




```julia
image(F)
```




    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <324, 558> with row [0 1 0 0]
    1//1 * <225, 104520> with row [0 0 1 0]




```julia
kernel(F)
```




    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]




```julia
is_injective(F), is_surjective(F)
```




    (true, false)



A map between mixed modules is of the data type `ModDedHom`. It is defined by a matrix of full rank over the corresponding number field. 


```julia
A = MixedMod(m1,m2);
B = MixedMod(m1,pseudo_matrix(matrix(k,4,4,[5 0 0 0 ; 0 6 0 0 ; 0 0 6 0 ; 0 0 0 8])));
Phi = ModDedHom(A,B,matrix(k,4,4,[5 0 0 0;0 6 0 0; 0 0 1 0; 0 0 0 1]))
```




    Module Map over Dedekind DomainReal quadratic field defined by x^2 - 10
    - - - -
    with domain:
    Mixed Module
    Mixed Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    - - -
    Torsion free part:
    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <3, 3> with row [1 0 0 0]
    1//1 * <1, 1> with row [0 1 0 0]
    - - -
    and torsion part:
    1-direct sum of quotient rings
    Base ring quotiented by the ideal <15, 165>
    - - - -
    and codomain
    Mixed Module
    Torsion Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    
    2-direct sum of quotient rings
    Base ring quotiented by the ideal <2, 2>
    Base ring quotiented by the ideal <120, 17160>
    - - - -
    defined by the matrix
    [5 0 0 0; 0 6 0 0; 0 0 1 0; 0 0 0 1]




```julia
image(Phi)
```




    Mixed Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    - - -
    Torsion free part:
    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <9, 48*sqrt(10) + 27> with row [1 0 0 0]
    1//1 * <1, 1> with row [0 1 0 0]
    - - -
    and torsion part:
    1-direct sum of quotient rings
    Base ring quotiented by the ideal <15, 165>




```julia
kernel(Phi)
```




    Mixed Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    - - -
    Torsion free part:
    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <12, 6*sqrt(10) + 42> with row [1 0 0 0]
    1//1 * <16, 8> with row [0 1 0 0]
    - - -
    and torsion part:
    1-direct sum of quotient rings
    Base ring quotiented by the ideal <30, 330>




```julia
Lambda = ModDedHom(A,B,matrix(k,4,4,[5 0 0 0;0 0 0 0; 0 0 1 0; 0 0 0 1]))
```

    Exception (fmpz_preinvn_init). Division by zero.





    Module Map over Dedekind DomainReal quadratic field defined by x^2 - 10
    - - - -
    with domain:
    Mixed Module
    Mixed Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    - - -
    Torsion free part:
    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <3, 3> with row [1 0 0 0]
    1//1 * <1, 1> with row [0 1 0 0]
    - - -
    and torsion part:
    1-direct sum of quotient rings
    Base ring quotiented by the ideal <15, 165>
    - - - -
    and codomain
    Mixed Module
    Torsion Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    
    2-direct sum of quotient rings
    Base ring quotiented by the ideal <2, 2>
    Base ring quotiented by the ideal <120, 17160>
    - - - -
    defined by the matrix
    [5 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1]




```julia
image(Lambda)
```

    Exception (fmpz_preinvn_init). Division by zero.
    Exception (fmpz_preinvn_init). Division by zero.





    Mixed Module over the base ring 
    Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    - - -
    Torsion free part:
    Pseudo-matrix over Maximal order of Real quadratic field defined by x^2 - 10 
    with basis nf_elem[1, sqrt(10)]
    1//1 * <9, 3> with row [1 0 0 0]
    1//1 * <1, 1> with row [0 1 0 0]
    - - -
    and torsion part:
    1-direct sum of quotient rings
    Base ring quotiented by the ideal <5, 5>




```julia
kernel(Lambda)
```

    Exception (fmpz_preinvn_init). Division by zero.



    Incompatible number of columns in vcat

    

    Stacktrace:

     [1] error(s::String)

       @ Base ./error.jl:33

     [2] vcat(a::AbstractAlgebra.Generic.MatSpaceElem{NfOrdElem}, b::AbstractAlgebra.Generic.MatSpaceElem{nf_elem})

       @ AbstractAlgebra ~/.julia/packages/AbstractAlgebra/ZmWFo/src/Matrix.jl:6328

     [3] pull_back(M::Hecke.PMat{nf_elem, Hecke.NfAbsOrdFracIdl{AnticNumberField, nf_elem}}, f::ModDedHom_tfree)

       @ Main ./In[12]:743

     [4] kernel(f::ModDedHom)

       @ Main ./In[12]:755

     [5] top-level scope

       @ In[13]:1

     [6] eval

       @ ./boot.jl:373 [inlined]

     [7] include_string(mapexpr::typeof(REPL.softscope), mod::Module, code::String, filename::String)

       @ Base ./loading.jl:1196
