#***************************************************************************
#                                   Reps
#       Copyright (C) 2005 Peter Webb 
#       Copyright (C) 2006 Peter Webb, Robert Hank, Bryan Simpkins 
#       Copyright (C) 2007 Peter Webb, Dan Christensen, Brad Froehle
#       Copyright (C) 2020 Peter Webb, Moriah Elkin
#       Copyright (C) 2022 Ferreira Juan David
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#  The overall structure of the reps package was designed and most of it
#  written by Peter Webb <webb@math.umn.edu>, who is also the maintainer. 
#  Contributions were made by
#  Dan Christensen, Roland Loetscher, Robert Hank, Bryan Simpkins,
#  Brad Froehle and others.
#***************************************************************************


InstallGlobalFunction( MakeIsoToPermGroup, function(rep)
    if IsBound(rep.permgroup) then
        return;
    fi;
    rep.isotopermgroup:=IsomorphismPermGroup(rep.group);
    rep.permgroup:=Image(rep.isotopermgroup,rep.group);
    return;
end );

InstallGlobalFunction( SafeNullspaceMat, function(mat, F)
    if mat=[] then return([]);
    fi;
    if mat[1]=[] then return(IdentityMat(Length(mat),F));
    fi;
    return(NullspaceMat(mat));
end );

#!
#! SafeBaseMat(M) returns a list a basis vectors for the space spanned
#! by the rows of M.  If M is an empty list, so M has no rows, or if
#! the elements of M are empty lists, so M has no columns, then the
#! empty list is returned.  
#!

InstallGlobalFunction( SafeBaseMat, function(M)
    if IsEmpty(M) or IsEmpty(M[1]) then
        return [];
    else
        return BaseMat(M);
    fi;
end );

#!
#! SafeIdentityMat(n, F) returns an nxn identity matrix over the field F.
#! When n = 0, returns an empty list.
#!

InstallGlobalFunction( SafeIdentityMat, function(n, F)
    if n=0 then
        return [];
    else
        return IdentityMat(n, F);
    fi;
end );


InstallGlobalFunction( SafeMatrixMult, function(A, B, n)
    if IsEmpty(A) then     # k = 0
        return [];
    elif n=0 then          # n = 0
        return List(A, row->[]);
    elif IsEmpty(B) then   # m = 0
        return NullMat(Length(A), n);
    else
        return A*B;
    fi;
end );

#!
#! SocleNullspaceMat(matrix, dimension, field)
#! returns a list of vectors forming a basis for the nullspace of the matrix,
#! and deals with the situation where the matrix may have no rows, etc.
#! It is used in SocleSeries(rep); SafeNullspaceMat works in other instances.
#! Written by Peter Webb July 2016.
#!
InstallGlobalFunction( SocleNullspaceMat, function(M,n,F)
    if n=0 or IsEmpty(M) then
        return SafeIdentityMat(n,F);
    else
        return NullspaceMat(M);
    fi;
end );
#!
#! InverseGenImages(rep) returns the list of matrices which are the inverses
#! of the matrices in rep.genimages. Two algorithms were tried: in
#! OldInverseGenImages, rather than invert the matrices, they
#! are raised to a power Order(g)-1 for each generator g. In fact, it appears
#! on testing this out that GAP is faster at doing x -> x^-1.
InstallGlobalFunction( InverseGenImages, function(rep)
    if IsBound(rep.inversegenimages) then
        return(rep.inversegenimages);
    fi;    
    rep.inversegenimages:=List(rep.genimages, x->x^-1);
    return(rep.inversegenimages);
end );

InstallGlobalFunction( OldInverseGenImages, function(rep)
    local orders, inverses, i;
    if IsBound(rep.inversegenimages) then
        return(rep.inversegenimages);
    fi;    
    orders:=List(GeneratorsOfGroup(rep.group),Order);
    inverses:=[];
    for i in [1..Length(orders)] do
        if orders[i]=infinity then
            inverses[i]:=rep.genimages[i]^-1;
        else
            inverses[i]:=rep.genimages[i]^(orders[i]-1);
        fi;
    od;
    rep.inversegenimages:=inverses;
    return(inverses);
end );


#!
#! RepToMatGroup(rep) returns the matrix group which is the image of the
#! representation.
#!


InstallGlobalFunction( RepToMatGroup, function(phi)
    return(Group(phi.genimages,IdentityMat(phi.dimension,phi.field)));
end );


InstallGlobalFunction( TrivialRep, function(group, field)
    local mats, onemat;
    onemat:=IdentityMat(1, field);
    mats:=List(GeneratorsOfGroup(group), g->onemat);
    return rec(
        group:=group,
        genimages:=mats,
        field:=field,
        dimension:=1,
        isRepresentation:=true,
        operations:=GroupRepOps
        );
end );


InstallGlobalFunction( TrivialRepN, function(group,field,n)
    local mats,mat;
    mat:=IdentityMat(n,field);
    mats:=List(GeneratorsOfGroup(group), g->mat);
    return rec(
        group:=group,
        genimages:=mats,
        field:=field,
        dimension:=1,
        isRepresentation:=true,
        operations:=GroupRepOps
        );
end );


InstallGlobalFunction( ZeroGroupRep, function(group,field)
    return rec(
     group:=group,
     genimages:=List(GeneratorsOfGroup(group), g->[]),
     field:=field,
     dimension:=0,
     isRepresentation:=true,
     operations:=GroupRepOps
     );
end );






#!
#!
#! DualRep(rep) 
#!
#! Updated Aug 2007 by Peter Webb using the function InverseGenImages.
#!


InstallGlobalFunction( DualRep , function(rep)
    if rep.dimension=0 then
        return(rep);
    fi;
    return rec(
     group:=rep.group,
     genimages:=List(InverseGenImages(rep), TransposedMat),
     field:=rep.field,
     dimension:=rep.dimension,
     isRepresentation:=true,
     operations:=GroupRepOps
     );
end );


#!
#! TensorProductMatrix(mat,mat) . . . Tensor product of two matrices
#!
#! The GAP command KroneckerProduct seems inexplicably slow and in tests on
#! two 60 x 60 matrices takes about twice as long as the following code.
#!


InstallGlobalFunction( TensorProductMatrix, function( A, B )
    local u, v, matrix;
    matrix := [ ];
    for u in A do
        for v in B do
            Add( matrix, Flat( List( u, x -> x * v ) ) );
        od;
    od;
    return( matrix );
end );

InstallGlobalFunction( TensorProductRep, function( rep1, rep2 )
    return rep1.operations.TensorProductRep( rep1, rep2 );
end );

InstallGlobalFunction( TensorProductGroupRep, function( rep1, rep2 )
    local mgens;
    if rep1.group <> rep2.group then
        Error("You must have two representations of the same group");
        return;
    fi;
    mgens := List( [1..Length( rep1.genimages )],
                    i -> TensorProductMatrix( rep1.genimages[i], rep2.genimages[i] ) );
    return Rep( rep1.group, mgens );
end );

InstallGlobalFunction( OldTensorProductRep, function(g,h)
    local mgens;
    if g.group <> h.group then
        Error("You must have two representations of the same group");
        return;
    fi;
    mgens := List( [1..Length( g.genimages ) ],
                   i -> KroneckerProduct( g.genimages[i], h.genimages[i] ) );
    return Rep(g.group,mgens);
end );


#!
#! TensorProductMorphism(M1,M2) Kronecker product of two morphisms
#!
#! Function introduced August 2007 by Dan Christensen.
#!


InstallGlobalFunction( TensorProductMorphism, function(M1,M2)
    return TensorProductMatrix(M1,M2);
end );

InstallGlobalFunction( OldTensorProductMorphism, function(M1,M2)
    return KroneckerProduct(M1,M2);
end );




InstallGlobalFunction( SubmoduleGroupRep, function(rep,v)
    local vs,base,newimages,g;
    vs:=VectorSpace(rep.field,v,[1..rep.dimension]*Zero(rep.field));
    base:=Basis(vs,v);
    newimages:=[];
    for g in rep.genimages do
        Add(newimages, List(base, b->Coefficients(base, b*g)));
    od;
    return rec(
        group:=rep.group,
        genimages:=newimages,
        field:=rep.field,
        dimension:=Length(base),
        isRepresentation:=true,
	operations:=GroupRepOps
	);
end );

#!
#! QuotientGroupRep is called by QuotientRep(rep, list of vecs) when rep is a
#! group representation.
#!
#!

InstallGlobalFunction( QuotientRep, function(rep,v)
    return rep.operations.QuotientRep(rep,v);
end );

InstallGlobalFunction( QuotientGroupRep, function(rep,v)
    local base,d,n,zero,onemat,i,positions,b,p,g,
        mat,newb,newimages,baseinverse, vs, tempbase;
    if Length(v)=0 then
        return rep;
    fi;
    base:=ShallowCopy(BaseMat(v));
    TriangulizeMat(base);
    d:=Length(base);
    n:=rep.dimension;
    if d=n then
        return ZeroGroupRep(rep.group,rep.field);
    fi;
    zero:=Zero(rep.field);
    onemat:=IdentityMat(n,rep.field);
    i:=1;
    positions:=[];
    for b in base do
        while b[i]=zero do
            Add(positions,i);
            i:=i+1;
        od;
        i:=i+1;
    od;
    Append(positions,[i..n]);
    for p in positions do
        Add(base, onemat[p]);
    od;
    baseinverse:=base^-1;
    newimages:=[];
    for g in rep.genimages do
        mat:=[];
        for p in positions do
            b:=g[p]*baseinverse;
            newb:=[];
            for i in [d+1..n] do
                Add(newb,b[i]);
            od;
            Add(mat,newb);
        od;
        Add(newimages, mat);
    od;
    return rec(
        group:=rep.group,
        genimages:=newimages,
        field:=rep.field,
        dimension:=n-d,
        isRepresentation:=true,
        operations:=GroupRepOps
        );
end );

#!
#! SectionRep(rep, list of vectors, list of vectors) returns the representation
#! on the quotient of the submodule spanned by the first list of vectors, by
#! the submodule spanned by the second list of vectors.
#! The two lists should be independent sets, and should be bases for
#! submodules of rep, the second submodule contained in the first.
#! This is not checked.
#! Written by Peter Webb July 2016.
#!
InstallGlobalFunction( SectionRep, function(rep,vectorsa,vectorsb)
    local repa, v, b, newvecsb;
    repa:=SubmoduleRep(rep,vectorsa);
    v:=VectorSpace(rep.field,vectorsa, [1..rep.dimension]*Zero(rep.field));
    b:=Basis(v,vectorsa);
    newvecsb:=List(vectorsb,x->Coefficients(b, x));
    return(QuotientRep(repa,newvecsb));
end );


#!
#!
#!PermutationMatrix(perm,d,field)  We specify d so that the permutation is
#! regarded as permuting [1..d]
#!
#!

InstallGlobalFunction( PermutationMatrix, function(perm,d, field)
    local m, p;
    m:=NullMat(d,d,field);
    for p in [1..d] do
        m[p][p^perm]:=One(field);
    od;
    return(m);
end );


#!
#!
#! PermToMatrixGroup( permgrp, field ) . . transforms a permutation group 
#! to a group of permutation matrices
#!
#!

#! PermToMatrixGroup , function( permgrp, field )
#! 
#!         local   matrix, g, p, d, mgens;
#! 
#!         d := LargestMovedPoint(permgrp);
#!         mgens:=[];
#!         for g in GeneratorsOfGroup(permgrp) do
#!                 matrix:=NullMat(d,d,field);
#!                 for p in [1..d] do
#!                         matrix[p][p^g]:=One(field);
#!                 od;
#!                 Add(mgens,matrix);
#!         od;
#!         return(Group(mgens));
#! end );

#!
#!
#! PermGroupToRep . . transforms a permutation group to a permutation
#!                                representation
#!

InstallGlobalFunction( PermGroupToRep , function( permgrp, field )
        local   matrix, g, p, d, mgens;
        d := LargestMovedPoint(permgrp);
        mgens:=[];
        for g in GeneratorsOfGroup(permgrp) do
            matrix:=NullMat(d,d,field);
            for p in [1..d] do
                matrix[p][p^g] := One(field);
            od;
            Add(mgens,matrix);
        od;
        return(Rep(permgrp, mgens));
end );

#!
#! Pretty print for matrices
#!

InstallGlobalFunction( DisplayMatrix, function( A )
    PrintArray( List (A, x -> List( x, Int ) ) );
end );

#!
#! PrintRep(rep) prints a representation nicely.
#!

InstallGlobalFunction( PrintRep, function( M )
    local g;
    if IsBound( M.name )  then
        Print( M.name );
    elif IsBound(M.longprint) and M.longprint then
        Print( "Representation( ", M.group, ", ", M.genimages, " )\n" );
    else
        Print( "Representation( ", M.group, ", Images \n");
        for g in M.genimages do
            DisplayMatrix(g);
        od;
        Print( " )\n" );
    fi;
end );


#!
#! CanRightRep(group,subgroup,list of elements, element)
#! returns the first element in the list which represents the same
#! coset as the element.
#! It does not verify the validity of its input.
#!


InstallGlobalFunction( CanRightRep, function(G,H,L,g)
   local h,k;
   h:=g^-1;
   for k in L do
     if k*h in H then return k;
     fi; 
   od;
end );


#!
#! RightCosetReps(group,subgroup) gives a list of representatives  of the 
#! right cosets of a group relative to a given subgroup.
#!


InstallGlobalFunction( RightCosetReps, function(G,H)
   return List(RightCosets(G,H), Representative);
end );

#! The function that gives the permutation representation of an element g
#! on a list L of right cosets of a subgroup has a built-in definition as
#! Permutation(g,L,OnRightCosets), but the problem is that L has to be a list 
#! whose elements are themselves lists, each containing the elements of
#! a Right Coset (Cosets are not lists).


#!
#! GLG(n,q) is a generalized version of GeneralLinearGroup
#! which does accept the case of dimension 1. 
#!


InstallGlobalFunction( GLG, function(n,q)
   if n=1 then
      return(Group(Z(q)*IdentityMat(1,GF(q))));
   else 
      return(GeneralLinearGroup(n,q));
   fi;
end );


#!
#! RepToGHBI(rep) turns a representation into a GroupHomomorphismByImages
#! with the corresponding data.
#!


InstallGlobalFunction( RepToGHBI, function(phi)
    local glg,h;
    glg := GLG(phi.dimension,Size(phi.field));
    h := GroupHomomorphismByImagesNC(
            phi.group,glg,
            GeneratorsOfGroup(phi.group),
            phi.genimages
        );
    #   h.isMapping:=true;
    #   h.isHomomorphism:=true;
    #   h.isGroupHomomorphism:=true;
    return h; 
end );



#!
#!
#! InducedRep(group, subgroup, rep of subgroup) computes the induction
#! of a rep from a subgroup to a group.
#!
#! In more detail, InducedRep(G, H, M) is M tensor_{kH} kG.
#! If m_1, ..., m_n is the standard basis for M and g_1, ..., g_r
#! are the coset representatives of the set {Hg} of right cosets,
#! then the basis chosen for the induced representation is
#!
#!   m_1 tensor g_1, m_2 tensor g_1, ..., m_n tensor g_1,
#!   m_1 tensor g_2, m_2 tensor g_2, ..., m_n tensor g_2,
#!   ...
#!   m_1 tensor g_r, m_2 tensor g_r, ..., m_n tensor g_r,
#!
#! in that order. 
#!


InstallGlobalFunction( InducedRep, function(G,H,phi)
    local L,n,F,R,Q,g,i,j,h,S,k,l,ghbiphi;
    ghbiphi:=RepToGHBI(phi);
    R:=[]; n:=phi.dimension; F:=phi.field;
    L:=RightCosetReps(G,H);
    for g in GeneratorsOfGroup(G) do
        Q:=NullMat(n*Length(L),n*Length(L),F);
        for i in [1..Length(L)] do
            j:=Position(L,CanRightRep(G,H,L,L[i]*g));
            h:=L[i]*g*(L[j]^-1);
            S:=ImagesRepresentative(ghbiphi,h);
            for k in [1..n] do
                for l in [1..n] do
                    Q[k+(i-1)*n][l+(j-1)*n]:=S[k][l];
                od;
            od;
        od;
    Add(R, Q);
    od;             
   return Rep(G,R);
end ); 


#!
#! InducedMorphism(group, subgroup, matrix) computes the induced
#! homomorphism of a map between representations of the subgroup.
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!


InstallGlobalFunction( InducedMorphism, function(G,H,S)
   local index;
   index := Index(G,H);
   return BlockMatrix(List([1..index], i->[i,i,S]), index, index);
   # Alternate implementation:
   # K := figure out the field
   # return KroneckerProduct(IdentityMat(index, K), S)
end );


#!
#! InducedInclusion(group, subgroup, rep of subgroup) computes the natural
#! homomorphism from rep to RestrictedRep(G, H, InducedRep(G, H, rep)).
#! With the basis conventions chosen here, this is just a matrix of the
#! form [ I | Z ], where I is an identity matrix and Z is a zero matrix.
#!
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!


InstallGlobalFunction( InducedInclusion, function(G,H,phi)
   local dim, index, M, F, i;
   F := phi.field;
   dim := phi.dimension;
   index := Index(G,H);
   M := NullMat(dim, dim*index, F);
   for i in [1..dim] do
      M[i][i] := One(F);
   od;
   return M;
end );


#!
#! InducedProjection(group, subgroup, rep of subgroup) computes the natural
#! projection from RestrictedRep(G, H, InducedRep(G, H, rep)) to rep.
#! With the basis conventions chosen here, this is just a matrix of the
#! form [ I ], where I is an identity matrix and Z is a zero matrix.
#!      [ Z ]
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!


InstallGlobalFunction( InducedProjection, function(G,H,phi)
   local dim, index, M, F, i;
   F := phi.field;
   dim := phi.dimension;
   index := Index(G,H);
   M := NullMat(dim*index, dim, F);
   for i in [1..dim] do
      M[i][i] := One(F);
   od;
   return M;
end );


#!
#! IsHom(rep1, rep2, mat) returns true if the matrix mat represents a
#! homomorphism from rep1 to rep2.  Quite useful for testing.
#!
#! Code written by Dan Christensen, University of Western Ontario,August 2007
#!


InstallGlobalFunction( IsHom , function(rep1, rep2, mat)
    local i;
    if rep1.group <> rep2.group then
        Error("You must have two representations of the same group");
        return;
    fi;
    for i in [1..Length(rep1.genimages)] do
        if rep1.genimages[i] * mat <> mat * rep2.genimages[i] then
            return false;
        fi;
    od;
    return true;
end );




InstallGlobalFunction( RegularRep, function(G, field)
    return PermutationRepOnCosets(G, Group(Identity(G)), field);
end );


InstallGlobalFunction( FreeRep, function(G, field, n)
    return TensorProductRep(TrivialRepN(G, field, n), RegularRep(G, field));
end );



#!
#! RelativeTrace(group,subgroup,rep) computes the matrix that corresponds to the mapping
#! tr_Q^P(v)=v(sum gi), where Q is a subgroup
#! of P, v is a vector in a P-representation phi, and the gi are a
#! complete list of representatives of right cosets of Q in P.
#!
#! The real use of this function is to apply the mapping to vectors which
#! are fixed by the subgroup, and then the result is fixed by the whole group.
#!
#! The algorithm is not very clever, and if the index of the subgroup is large
#! it would be better to construct a chain of subgroups between the two groups
#! and compute the relative trace as the product of the relative traces between
#! pairs of groups in the chain. Such an algorithm is used in the command
#! NormRep which returns the relative trace from the identity subgroup.
#!


InstallGlobalFunction( RelativeTrace, function(P,Q,rep)
   local M,l,g,ghbi;
   M:=NullMat(rep.dimension,rep.dimension,rep.field);
   l:=RightCosetReps(P,Q);
   ghbi:=RepToGHBI(rep );
   for g in l do
      M:=M+Image(ghbi,g);
   od;
   return M;
end );

#!
#! OldHomBasis(rep1,rep2) . . returns a basis for the space of 
#! module homomorphisms A -> B. The elements of this basis are matrices.
#!
#! This function is not as fast as the code for HomBasis written by 
#! Dan Christensen, which has a more straightforward setup matrix whose
#! nullspace we find.
#!


InstallGlobalFunction( OldHomBasis, function(g,h)
local  f, basis, i, j, m, u;
f:=FixedPoints@(TensorProductRep(DualRep(g),h));
basis:=[];
u:=[];
for m in f do
    for i in [1..g.dimension] do
    u[i]:=[];
        for j in [1..h.dimension] do
            u[i][j]:=m[h.dimension*(i-1)+j];
        od;
    od;
    Add(basis,ShallowCopy(u));
od;
return(basis);
end );


#!
#! HomBasis(M, N) returns a basis for the space of 
#! kG-module homomorphisms M -> N. The elements of this basis are matrices.
#!
#! The algorithm used finds matrices X so that Xg-gX=0 for all group generators
#! g. This sets up a linear algebra problem more efficiently than solving
#! gXg^-1=X, and was observed independently by John Pliam (1993),
#! Michael Smith (1993) and Dan Christensen (2007).
#! This code written by Dan Christensen, August 2007.
#!
#! GroupHomBasis is called by HomBasis(M, N) when M is a
#! group representation.
#!


InstallGlobalFunction( HomBasis, function(M,N)
    return M.operations.HomBasis(M,N);
end );

InstallGlobalFunction( GroupHomBasis, function(M,N)
    local dimM, dimN, r, i, j, k, l, v, f, basis, m, u;
    dimM := M.dimension;
    dimN := N.dimension;
    r := Length(M.genimages);
    v:=NullMat(dimM*dimN, dimM*dimN*r, M.field);
    for i in [1..dimM] do
        for k in [1..dimN] do
            for l in [1..r] do
                for j in [1..dimM] do
                    v[(j-1)*dimN+k][(l-1)*dimM*dimN+(i-1)*dimN+k] :=
                    v[(j-1)*dimN+k][(l-1)*dimM*dimN+(i-1)*dimN+k] + M.genimages[l][i][j];
                od;
                for j in [1..dimN] do
                    v[(i-1)*dimN+j][(l-1)*dimM*dimN+(i-1)*dimN+k] :=
                    v[(i-1)*dimN+j][(l-1)*dimM*dimN+(i-1)*dimN+k] - N.genimages[l][j][k];
                od;
            od;
        od;
    od;
    f := NullspaceMat(v);
    basis:=[];
    for m in f do
        u:=[];
        for i in [1..M.dimension] do
            u[i]:=[];
            for j in [1..N.dimension] do
                u[i][j]:=m[N.dimension*(i-1)+j];
            od;
        od;
        Add(basis, u);
    od;
    return(basis);
end );


#!
#! DimHom(rep1,rep2) . . returns the dimension of the space of 
#! module homomorphisms rep1 -> rep2.
#!


InstallGlobalFunction( DimHom, function(g,h)
    return Length(HomBasis(g,h));
end );


#!
#! SumOfImages(rep,rep) . . returns a basis for the sum of images of all 
#! module homomorphisms A -> B
#!
#! GroupSumOfImages is called by SumOfImages(rep,rep) when rep is a
#! group representation.
#!
#! Updated by Peter Webb Sep 11, 2007
#!


InstallGlobalFunction( SumOfImages, function(M,N)
    return M.operations.SumOfImages(M,N);
end );

InstallGlobalFunction( GroupSumOfImages, function(g,h)
    return SafeBaseMat(Concatenation(HomBasis(g,h)));
end );



#!
#! DecomposeSubmodule(repn, basis) . . probably returns a list of two bases 
#! for summands of 
#! the submodule of repn spanned by basis, if the module is decomposable, and
#! if the module is indecomposable it returns a list whose only element is the
#! given basis.
#!
#! GroupDecomposeSubmodule is called by DecomposeSubmodule(repn, basis) when repn
#! is a group representation.
#!
#! This code was written by Bryan Simpkins and Robert Hank, University of
#! Minnesota, April 2006. A related approach was taken by Michael Smith in
#! code written in 1993.
#! It was tidied up by Dan Christensen, University of Western Ontario, Aug 2007.
#! In July 2016 Craig Corsi, University of Minnesota, made a change so that
#! decompositions over fields larger than the field of definition of the 
#! representation are found.
#!


InstallGlobalFunction( DecomposeSubmodule, function(M,basisstructure)
    return M.operations.DecomposeSubmodule(M,basisstructure);
end );

InstallGlobalFunction( GroupDecomposeSubmodule, function(repn, basis)
	local newrep, initlist, a, kernel, image, b, d, n, z;
	if Length(basis) <= 1 then return [basis]; fi;
	newrep := SubmoduleRep(repn, basis);
	z:=PrimitiveElement(repn.field);
	initlist := HomBasis(newrep,newrep);
	Add(initlist, 0 * initlist[1]);
	b := Length(initlist);
	while b > 0 do;
		d := b - 1;
		while d > 0 do;
			a := z*initlist[b] + initlist[d];
			n := 1;
			while n < Length(basis) do;
				a:= a*a;
				n:= 2*n;
			od;
			kernel := NullspaceMat(a);
			if not(Length(kernel) = 0 or Length(kernel) = Length(basis)) then
				image  := BaseMat(a);
				return [kernel * basis, image * basis];
			fi;
			d := d - 1;
		od;
		b := b - 1;
	od;
	return [basis];
end );


#!
#! DecomposeOnce(rep) . . probably returns a list of two bases for summands
#! of the representation rep if it is decomposable, and returns a list whose
#! only element is the standard basis if it is indecomposable. 
#!
#! This code was written by Bryan Simpkins and Robert Hank, University of
#! Minnesota, April 2006
#!


InstallGlobalFunction( DecomposeOnce, function(rep)
	return DecomposeSubmodule(rep, IdentityMat(rep.dimension, rep.field));
end );


#!
#! Decompose(rep) . . returns a list of bases of direct summands of rep which
#! are probably indecomposable.
#!
#! DecomposeGroupRep is called by Decompose(rep) when rep
#! is a group representation.
#!
#! This code was written by Bryan Simpkins and Robert Hank, University of
#! Minnesota, April 2006.
#! It was modified by Dan Christensen and Peter Webb August 2007 so that indecomposable
#! representations are only checked once (before they were checked twice),
#! so as to use space more economically, and so that the result is stored
#! in rep.summands.
#!


InstallGlobalFunction( Decompose, function(rep)
    return rep.operations.Decompose(rep);
end );

InstallGlobalFunction( DecomposeGroupRep, function(rep)
	local summands, result, q;
	if IsBound(rep.summands) then return(rep.summands);
	fi;
	summands := [IdentityMat(rep.dimension,rep.field)];
	q := 1;
	# We maintain the following invariants:
	# - summands is a list of lists of vectors; the union of these
	#   lists forms a basis for rep.field^(rep.dimension).
        # - the summands at positions < q appear to be indecomposable;
	#   those at positions >= q haven't been investigated completely.
	while IsBound(summands[q]) do;
		result := DecomposeSubmodule(rep, summands[q]);
		if Length(result) = 2 then
			summands[q] := result[1];
			Add(summands, result[2]);
		else
			q := q + 1;
		fi;
	od;
	rep.summands:=summands;
	return summands;
end );


#!
#! MatsOfCosetReps(rep) . . returns a list
#!
#! [rep of a proper subgroup, 
#! list of matrices of right coset representatives of the subgroup,
#! list of right coset representatives of the subgroup].
#!
#! The coset representatives are a Schreier
#! transversal for the stabilizer of a point, and the stabilizer subgroup has Schreier
#! generators for the stabilizer. At the same time matrices which
#! represent the elements which arise are stored.
#!
#! This function is called by NormRep, MatrixOfElement and RestrictedRep.
#! 


InstallGlobalFunction( MatsOfCosetReps, function(rep)
    local s, orbit, positions, mats, inversemats, perms, gens, q, x, j, ximage, table,
    subgroupgens, subgroupgenimages;
    s:=SmallestMovedPoint(rep.group);
    orbit:=[s];
    positions:=[];
    mats:=[];
    mats[s]:=IdentityMat(rep.dimension,rep.field);
    inversemats:=[];
    inversemats[s]:=IdentityMat(rep.dimension,rep.field);
    perms:=[];
    perms[s]:=();
    gens:=GeneratorsOfGroup(rep.group);
    InverseGenImages(rep);
    table:=[];
    for j in [1..Length(gens)] do
        table[j]:=[];
    od;
    q:=1;
    while IsBound(orbit[q]) do
        x:=orbit[q];
        for j in [1..Length(gens)] do
            ximage:=x^gens[j];
            if not ximage in orbit then
                Add(orbit,ximage);
                positions[ximage]:=j;
                perms[ximage]:=perms[x]*gens[j];
                mats[ximage]:=mats[x]*rep.genimages[j];
                inversemats[ximage]:=rep.inversegenimages[j]*inversemats[x];
                else table[j][x]:=perms[x]*gens[j]*perms[ximage]^-1;
            fi;
        od;
        q:=q+1;
    od;
    subgroupgens:=[];
    subgroupgenimages:=[];
    for j in [1..Length(gens)] do
        for x in orbit do
            if IsBound(table[j][x]) and
            Size(Group(subgroupgens,()))<
            Size(Group(Concatenation(subgroupgens,[table[j][x]])))
                then Add(subgroupgens,table[j][x]);
                Add(subgroupgenimages,mats[x]*rep.genimages[j]*inversemats[x^gens[j]]);
            fi;
        od;
    od;
    return ([rec(
        group:=Group(subgroupgens,()),
        genimages:=subgroupgenimages,
        field:=rep.field,
        dimension:=rep.dimension,
        isRepresentation:=true,
        operations:=GroupRepOps
        ),
        mats,perms]);
end );


#!
#! NormRep(rep) . . returns the matrix which represents the sum of the group
#!                  elements.
#!
#! NormRep calls MatsOfCosetReps recursively until the subgroup is the identity group, and
#! multiplies the sums of the matrices of the coset representatives. 
#! 


InstallGlobalFunction( NormRep, function(rep)
    local output, newrep, temp, a, sum;
    output:=IdentityMat(rep.dimension,rep.field);
    MakeIsoToPermGroup(rep);
    newrep:=rec(
        group:=rep.permgroup,
        genimages:=rep.genimages,
        field:=rep.field,
        dimension:=rep.dimension,
        isRepresentation:=true,
        operations:=GroupRepOps
    );
    while Size(newrep.group)>1 do
        temp:=MatsOfCosetReps(newrep);
        newrep:=temp[1];
        sum:=NullMat(rep.dimension,rep.dimension,rep.field);
        for a in temp[2] do
            sum:=sum+a;
        od;
        output:=sum*output;
    od;
    return output;
end );


InstallGlobalFunction( OldNormRep, function(rep)
        local p, currentgroup, stablist, subgroup, replist, norm;
        currentgroup := rep.group;
        p := RepToGHBI(rep);
        norm := IdentityMat(rep.dimension,rep.field);
        while not IsTrivial(currentgroup) do;
             stablist := StabChain(currentgroup);
             subgroup := Subgroup(currentgroup,stablist.stabilizer.generators);
             replist := RightCosetReps(currentgroup, subgroup);
             replist := List( replist, x -> ImagesRepresentative(p, x));
             norm := Sum(replist) * norm;
             currentgroup := subgroup;
        od;
        return norm;
end );


#!
#! MatrixOfElement(rep,group element) returns the matrix which represents
#! the group element.
#!


InstallGlobalFunction( MatrixOfElement, function(rep,g)
local newrep, newg, mat, s, temp;
    MakeIsoToPermGroup(rep);
    newrep:=rec(
        group:=rep.permgroup,
        genimages:=rep.genimages,
        field:=rep.field,
        dimension:=rep.dimension,
        isRepresentation:=true,
        operations:=GroupRepOps
    );
    newg:=Image(rep.isotopermgroup,g);
    mat:=IdentityMat(rep.dimension,rep.field);
    while Size(newrep.group)>1 do
        s:=SmallestMovedPoint(newrep.group);
        temp:=MatsOfCosetReps(newrep);
        newrep:=temp[1];
        mat:=temp[2][s^newg]*mat;
        newg:=newg*temp[3][s^newg]^-1;
    od;
    return(mat);
end );


#!
#! MatricesOfElements(rep,list of group elements) returns the list of matrices 
#! which represent the group elements.
#!


InstallGlobalFunction( MatricesOfElements, function(rep,l)
local newrep, newl, matlist, s, temp, x, n;
    MakeIsoToPermGroup(rep);
    newrep:=rec(
     group:=rep.permgroup,
     genimages:=rep.genimages,
     field:=rep.field,
     dimension:=rep.dimension,
     isRepresentation:=true,
     operations:=GroupRepOps
     );
    n:=Length(l);
    newl:=List(l, x->Image(rep.isotopermgroup,x));
    matlist:=List(l,x->IdentityMat(rep.dimension,rep.field));
    while Size(newrep.group)>1 do
        s:=SmallestMovedPoint(newrep.group);
        temp:=MatsOfCosetReps(newrep);
        newrep:=temp[1];
        matlist:=List([1..n], x->temp[2][s^newl[x]]*matlist[x]);
        newl:=List(newl,x->x*temp[3][s^x]^-1);
    od;
    return(matlist);
end );


#!
#! RestrictedRep(group, subgroup, rep of group) computes the restriction
#! of a rep from a group to a subgroup.
#!
#! This code has been rewritten September 2007 by Peter Webb.
#! The algorithm calls the function MatricesOfElements which finds the
#! matrices representing the generators of the subgroup by working down a
#! stabilizer chain for the group and at each stage computing matrices
#! which represent the coset representatives of the stabilizer subgroups
#! and also the generators of the stabilizers. This approach is more 
#! economical that expressing each subgroup generator as a word in the 
#! generators of the big group and then evaluating that word on matrices, because
#! in such words there are repeated subwords which get evaluated again and
#! again.
#!
#! The previous implementation of this function appears as OldRestrictedRep.
#!
#!

InstallGlobalFunction( RestrictedRep, function(G,H,rep)
    local R;
    R:=MatricesOfElements(rep,GeneratorsOfGroup(H));
    return rec(
        group:=H,
        genimages:=R,
        field:=rep.field,
        dimension:=rep.dimension,
        isRepresentation:=true,
        operations:=GroupRepOps
        );
end );



InstallGlobalFunction( OldRestrictedRep, function(G,H,phi)
   local ghbiphi,R;
   ghbiphi:=RepToGHBI(phi);
   R:=List(GeneratorsOfGroup(H), g->ImagesRepresentative(ghbiphi, g));
    return rec(
        group:=H,
        genimages:=R,
        field:=phi.field,
        dimension:=phi.dimension,
        isRepresentation:=true,
        operations:=GroupRepOps
        );
end );



#!
#! SymmetricPowerRep(rep, n) computes the action of rep on the n-th symmetric power
#! of the vector space V associated with rep.
#! This routine was written by Brad Froehle at the University of Minnesota, January 2007.
#!


InstallGlobalFunction( SymmetricPowerRep , function(rep, n)
	local Sn;

	#! Define a function Sn which can be called repeatedly for each matrix corresponding to a generator of the group.
	Sn := function(A, n)
		local MyProduct, DualList, dimV, expEnd, dimOut, output, waysToSum, dual, i, j, l;

		#! MyProduct(A, x, y)
		#! 
		#! Assumptions: A is a square n-by-n matrix, x and y are lists of the same length
		#! whose entries are elements of [1..n]
		#!
		#! Output: the product A[x[1]][y[1]] * A[x[2]][y[2]] * A[x[3]][y[3]] * ...
		#! 
		MyProduct := function(A, x, y)
			local output, i;
			output := 1;
			for i in [1..Length(x)]
				do
				output := output * A[x[i]][y[i]];
			od;
			return output;
		end;
		
		#! DualList(A) 
		#!
		#! Assumptions: A is a list.
		#!
		#! Output: A list in which 1 is repeated A[1] times, 2 is repeated A[2] times, etc...
		#!
		DualList := function(A)
			local i, output;
			output := [];
			for i in [1..Length(A)]
				do
				Append(output, List([1..A[i]], x->i));
			od;
			return output;
		end;

		#! Define some variables which are used repeatedly:
		#!  * dimV is the dimension of the initial representation.
		#!  * expEnd is a list whose elements describe all possible monomials of total degree dimV.
		#!  * dimOut is the dimension of the output representation, i.e. the number of monomials in expEnd.
		#!
		dimV := Size(A);
		expEnd := Reversed(OrderedPartitions(dimV+n,dimV)-1);
		dimOut := Size(expEnd);
		
		#! Initialize output to a dimOut x dimOut matrix of all zeros
		#!
		output := [];
		for i in [1..dimOut]
			do
			output[i] := List([1..dimOut],x->0);
		od;

		#! Iterate, calculating each entry of output individually.
		#!
		for i in [1..dimOut]
			do
			# Because this next call (PermutationsList) might be slow, we calculate it as few times as possible.
			waysToSum := PermutationsList(DualList(expEnd[i]));
			for j in [1..dimOut]
				do
				# Calculate the ji-th entry here
				dual := DualList(expEnd[j]);
				for l in waysToSum
					do
					output[j][i] := output[j][i] + MyProduct(A, dual, l);
				od;
			od;
		od;
		return output;
	end;
	#! End function Sn.
	
	return Rep(rep.group, List(rep.genimages, x-> Sn(x,n)));
end );


#!
#! ProjectiveHomBasis(rep1,rep2) returns a list of matrices which form a basis for
#! the space of module homomorphisms from rep1 to rep2 which factor through
#! a projective module. The algorithm computes the image of the norm map
#! applied to the representation Dual(rep1) tensor rep2.
#!


InstallGlobalFunction( ProjectiveHomBasis, function( rep1, rep2 )
    local f, basis, i, j, m, u;
    f:=SafeBaseMat( NormRep( TensorProductRep( DualRep( rep1 ), rep2 ) ) );
    basis:=[];
    for m in f do
        u:=[];
        for i in [1..rep1.dimension] do
        u[i]:=[];
            for j in [1..rep2.dimension] do
                u[i][j]:=m[rep2.dimension*(i-1)+j];
            od;
        od;
        Add(basis,u);
    od;
    return(basis);
end );

InstallGlobalFunction( OldProjectiveHomBasis, function(rep1,rep2)
    local  f, basis, i, j, m, u;
    f:=SafeBaseMat(OldNormRep(TensorProductRep(DualRep(rep1),rep2)));
    basis:=[];
    u:=[];
    for m in f do
        for i in [1..rep1.dimension] do
        u[i]:=[];
            for j in [1..rep2.dimension] do
                u[i][j]:=m[rep2.dimension*(i-1)+j];
            od;
        od;
        Add(basis,ShallowCopy(u));
    od;
    return(basis);
end );



#!
#! IsProjectiveMorphism(rep1, rep2, f) returns whether or not
#! the morphism f: rep1 --> rep2 factors through a projective.
#! f is assumed to be a kG-module homomorphism.
#!


InstallGlobalFunction( IsProjectiveMorphism, function(M,N,f)
    local mat, n;
    mat:=ProjectiveHomBasis(M,N);
    n:=Length(mat);
    Add(mat,f);
    mat:=List(mat,Flat);
    return RankMat(mat)=n;
end );



#!
#! IsProjectiveRep(rep) returns true if the representation is a projective
#! module, and false otherwise. The algorithm restricts the representation
#! to a Sylow p-subgroup and tests whether |G| times the rank of the norm
#! map equals the dimension of the representation.
#!


InstallGlobalFunction( IsProjectiveRep, function(rep)
    local s, resrep, n;
    if Characteristic(rep.field)=0 then return(true);
    fi;
    s:=SylowSubgroup(rep.group, Characteristic(rep.field));
#    s:=Subgroup(rep.group,SmallGeneratingSet(s));  This line is for permutation groups
    resrep:=RestrictedRep(rep.group, s, rep);
    n:=NormRep(resrep);
    if Rank(n)*Size(s) = rep.dimension then
        return(true);
    else return(false);
    fi;
end );


#!
#! ProjectiveFreeSummand(rep) returns a summand of rep which is probably 
#! projective free.  The complementary summand is guaranteed to be 
#! projective, so the returned rep is stably isomorphic to the original rep.
#!


InstallGlobalFunction( ProjectiveFreeSummand, function(M)
    local comps, comp, pfbasis;
    comps:=Decompose(M);
    pfbasis:=[];
    for comp in comps do
        if not IsProjectiveRep(SubmoduleRep(M, comp)) then
            Append(pfbasis, comp);
        fi;
    od;
    return SubmoduleRep(M, pfbasis);
end );


#!
#! ProjectiveDecomposition(rep) returns a list of two bases, for a submodule
#! which is projective, and for a submodule with probably no non-zero projective summands,
#! whose direct sum is the whole representation. 
#!


InstallGlobalFunction( ProjectiveDecomposition, function(M)
    local comps, comp, pfbasis, pbasis;
    comps:=Decompose(M);
    pfbasis:=[];
    pbasis:=[];
    for comp in comps do
        if not IsProjectiveRep(SubmoduleRep(M, comp)) then
            Append(pfbasis, comp);
        else Append(pbasis, comp);
        fi;
    od;
    return ([pbasis,pfbasis]);
end );

#!
#! Interface to the Meataxe.
#!
#! In the following routines several of the meataxe commands available in GAP
#! are converted so that they take representations as arguments and return
#! representations, where the meataxe implementation would have returned
#! meataxe modules. Clearly more of the meataxe commands could be implemented
#! than is done below.
#!

#!
#! MeataxeModuleToRep(rep,meataxemodule) converts a meataxe module to a
#! representation. Because the meataxe module does not store the group
#! being represented, this is obtained from a representation called rep
#! whose only property is that rep.group is the required group.
#!


InstallGlobalFunction( MeataxeModuleToRep, function(rep,meataxemodule)
    return(rec(group:=rep.group,
           genimages:=MTX.Generators(meataxemodule),
           field:=MTX.Field(meataxemodule),
           dimension:=MTX.Dimension(meataxemodule),
           isRepresentation:=true,
           operations:=GroupRepOps)
        );
end );


#!
#! RepToMeataxeModule(rep) converts a representation to a meataxe module.
#!


InstallGlobalFunction( RepToMeataxeModule, function(rep)
    if rep.genimages=[] then
        return(GModuleByMats([], rep.dimension, rep.field ));
    elif rep.dimension=0 then
        return(GModuleByMats([], rep.dimension, rep.field ));
    else
        return(GModuleByMats(rep.genimages, rep.field ));
    fi;
end );


#!
#! ProperSubmodule(rep) calls the meataxe command MTX.ProperSubmoduleBasis
#! and returns a basis for a proper submodule, or [] if there is none.
#!


InstallGlobalFunction( ProperSubmodule, function(rep)
    local basis;
    if rep.dimension=0 then return([]); fi;
    basis:=(MTX.ProperSubmoduleBasis(RepToMeataxeModule(rep)));
    if basis=fail then return([]); fi;
    return(basis);
end );

InstallGlobalFunction( IsIrreducibleRep, function(rep)
    if rep.dimension=0 then return(false); fi;
    return(MTX.IsIrreducible(RepToMeataxeModule(rep)));
end );

InstallGlobalFunction( IsAbsolutelyIrreducibleRep, function(rep)
    if rep.dimension=0 then return(false); fi;
    return(MTX.IsAbsolutelyIrreducible(RepToMeataxeModule(rep)));
end );

InstallGlobalFunction( BasesCompositionSeriesRep, function(rep)
    if rep.dimension=0 then return([[]]); fi;
    return(MTX.BasesCompositionSeries((RepToMeataxeModule(rep))));
end );

InstallGlobalFunction( CompositionFactorsRep, function(rep)
    local modules;
    if rep.dimension=0 then return([]); fi;
    modules:=MTX.CompositionFactors(RepToMeataxeModule(rep));
    return(List(modules,x->MeataxeModuleToRep(rep,x)));    
end );

InstallGlobalFunction( RadicalRep, function(rep)
    if rep.dimension=0 then return([]); fi;
    return(MTX.BasisRadical(RepToMeataxeModule(rep)));
end );

InstallGlobalFunction( SocleRep, function(rep)
    if rep.dimension=0 then return([]); fi;
    return(MTX.BasisSocle(RepToMeataxeModule(rep)));
end );



#!
#!
#! FixedPointRep(representation, subgroup of rep.group) returns the 
#! representation of the normalizer of the subgroup on the fixed points
#! of the subgroup.  Written by Peter Webb June 2016.
#!
#!

InstallGlobalFunction( FixedPointRep, function(rep,gp)
    local n, resn, resgp;
    n:=Normalizer(rep.group,gp);
    resn:=RestrictedRep(rep.group,n,rep);
    resgp:=RestrictedRep(rep.group,gp,rep);
    return(SubmoduleRep(resn,FixedPoints@(resgp)));
end );

#!
#!
#! BrauerRep(representation, p-subgroup of rep.group) returns the 
#! representation of the normalizer of the subgroup on the fixed points
#! of the subgroup modulo the image of traces from proper subgroups of 
#! the p-group.  Written by Peter Webb June 2016.
#!
#!

InstallGlobalFunction( BrauerRep, function(rep,p)
    local n, resn, resp, fp, fprep, maxsub, traceimages, h, M, resh, v, b;
    n:=Normalizer(rep.group,p);
    resn:=RestrictedRep(rep.group,n,rep);
    resp:=RestrictedRep(rep.group,p,rep);
    fp:=FixedPoints@(resp);
    fprep:=SubmoduleRep(resn,fp);
    maxsub:=MaximalSubgroups(p);
    traceimages:=[];
    for h in maxsub do
        M:=RelativeTrace(p,h,resp);
        resh:=RestrictedRep(resp.group,h,resp);
        Append(traceimages,SafeMatrixMult(FixedPoints@(resh),M,rep.dimension));
    od;
    traceimages:=SafeBaseMat(traceimages);
    v:=VectorSpace(rep.field,fp);
    b:=Basis(v,fp);
    traceimages:=List(traceimages,x->Coefficients(b, x));
    return QuotientRep(fprep, traceimages);
end );

#!
#!
#! RadicalSeries(rep) returns a list with two entries [bases, reps] where 
#! bases is a list of bases of the successive powers of the radical
#! rep = rad^0(rep), rad^1(rep), ...
#! in descending order. The last term in the list is the empty list.
#! reps is a list of the representations on the radical quotients
#! rep/rad^1(rep), rad^1(rep)/rad^2(rep). ...
#! all of which are semisimple. The last term is the last nonzero 
#! representation, and so the list is one shorter than bases.
#! Written by Peter Webb July 2016.
#!
#!

InstallGlobalFunction( RadicalSeries, function(rep)
    local rad, i, submod, submodrad, layer;
    rad:=[];
    layer:=[];
    rad[1]:=SafeIdentityMat(rep.dimension,rep.field);
    i:=1;
    while Length(rad[i])>0 do
        i:=i+1;
        submod:=SubmoduleRep(rep,rad[i-1]);
        submodrad:=RadicalRep(submod);
        rad[i]:=SafeMatrixMult(submodrad,rad[i-1],submod.dimension);
        layer[i-1]:=QuotientRep(submod,submodrad);
    od;
    return([rad,layer]);
end );

#!
#!
#! SocleSeries(rep) returns a list with two entries [bases, reps] where 
#! bases is a list of bases of the higher socles
#! rep = soc^t(rep), soc^(t-1)(rep), ...
#! in DESCENDING order. The last term in the list is the empty list.
#! reps is a list of the representations on the socle quotients
#! rep/soc^(t-1)(rep), soc^(t-1)(rep)/soc^(t-2)(rep). ...
#! all of which are semisimple. The last term is the last nonzero 
#! representation, and so the list is one shorter than bases.
#! Written by Peter Webb July 2016.
#!
#!
InstallGlobalFunction( SocleSeries, function(rep)
    local radseries, n, temp, socseries, i;
    radseries:=RadicalSeries(DualRep(rep));
    socseries:=[];
    socseries[1]:=[];
    socseries[2]:=[];
    temp:=[];
    temp[1]:=List(radseries[1],y->SocleNullspaceMat(TransposedMat(y),rep.dimension,rep.field));
    temp[2]:=List(radseries[2],DualRep);
    n:=Length(radseries[2]);
    for i in [1..n+1] do
        socseries[1][i]:=temp[1][n+2-i];
    od;
    for i in [1..n] do
        socseries[2][i]:=temp[2][n+1-i];
    od;
    return(socseries);
end );

#!
#! ButterflyFactors(rep, descending filtration, descending filtration) 
#! returns a matrix whose entries are the representations that appear as
#! sections in Zassenhaus' Butterfly Lemma.
#! Each descending filtration is a list of bases of submodules of rep, forming
#! a descending chain. The representations in the output are the factors in
#! a common refinement of the two filtrations, fand their position in the 
#! refinement is indicated by their position in the matrix.
#! Written by Peter Webb July 2016
#!
InstallGlobalFunction( ButterflyFactors, function(rep,L1,L2)
    local factors, i, j, vecsa, vecsb;
    factors:=[];
    for i in [1..Length(L1)-1] do
        factors[i]:=[];
        for j in [1..Length(L2)-1] do
            vecsa:=SumIntersectionMat(L1[i],L2[j])[2];
            vecsb:=SumIntersectionMat(SumIntersectionMat(L1[i+1],L2[j])[2],
            SumIntersectionMat(L1[i],L2[j+1])[2])[1];
            factors[i][j]:=SectionRep(rep,vecsa,vecsb);
        od;
    od;
    return(factors);
end );

InstallGlobalFunction( DirectSumRep, function(rep1,rep2)
    return rep1.operations.DirectSumRep(rep1,rep2);
end );


#!
#! DirectSumGroupRep(rep1,rep2) returns the representation that is the direct sum of 
#! representations rep1 and rep2 for the same group. Written by Peter Webb February 2020.
#!
#! DirectSumGroupRep is called by DirectSumRep(rep1,rep2) when rep1 and rep2 are
#! group representations.
#!
#!


InstallGlobalFunction( DirectSumGroupRep, function( rep1, rep2 )
    local i, gimages, x, zerovec1, zerovec2, padded1, padded2;
    if rep1.field<>rep2.field or rep1.group<>rep2.group then
        Error("Representations incompatible.");
    fi;
    zerovec1:=Zero(rep1.field)*[1..rep1.dimension];
    zerovec2:=Zero(rep2.field)*[1..rep2.dimension];
    gimages:=[];
    for i in [1..Length(rep1.genimages)] do
    padded1:=List(rep1.genimages[i],x->Concatenation(x,zerovec2));
    padded2:=List(rep2.genimages[i],x->Concatenation(zerovec1,x));
    gimages[i]:=Concatenation(padded1,padded2);
    od;
    return rec(
               group:=rep1.group,
               genimages:=gimages,
               field:=rep1.field,
               dimension:=rep1.dimension+rep2.dimension,
               isRepresentation:=true,
               operations:=GroupRepOps
            );
end );


#!
#! ProjectiveFreeCore(rep) returns the representation that is the largest direct summand of 
#! rep with not projective summand. It is only guaranteed to work when the group is a p-group
#! in characteristic p. In other cases it may give a correct answer, and if it does not then an error
#! message is returned. Written by Peter Webb February 2020. There is another function 
#! ProjectiveFreeSummand which invokes Decompose, and because of this will be more limited
#!


InstallGlobalFunction( ProjectiveFreeCore, function(rep)
    local n, resrep, vectors, sub;
    resrep:=RestrictedRep(rep.group, SylowSubgroup(rep.group,Characteristic(rep.field)),rep);
    n:=NormRep(resrep);
    if Length(n)=0 then
        return(rep);
    fi;
    vectors:=BaseSteinitzVectors(IdentityMat(rep.dimension,rep.field),NullspaceMat(n));
    sub:=Spin(rep,vectors.factorspace);
    if IsProjectiveRep(SubmoduleRep(rep,sub))=false then
        Error("The routine tried to factor out a non-projective module. It only works for p-groups.");
        return;
    fi;
    return(QuotientRep(rep,sub));
end );

#!
#! ProjectiveSummand(rep) returns a basis for a maximal projective direct summand of
#! rep. It is only guaranteed to work when the group is a p-group in characteristic p.
#! In other cases it may give a correct answer, and if it does not then an error
#! message is returned. It does not use Decompose. Written by Peter Webb February 2020.
#!

InstallGlobalFunction( ProjectiveSummand, function(rep)
    local n, resrep, vectors, sub;
    resrep:=RestrictedRep(rep.group, SylowSubgroup(rep.group,Characteristic(rep.field)),rep);
    n:=NormRep(resrep);
    if Length(n)=0 then
        return(rep);
    fi;
    vectors:=BaseSteinitzVectors(IdentityMat(rep.dimension,rep.field),NullspaceMat(n));
    sub:=Spin(rep,vectors.factorspace);
    if IsProjectiveRep(SubmoduleRep(rep,sub))=false then
        Error("The routine produced a basis called sub, for a submodule that is too big and not projective. It is only guaranteed to work for p-groups.");
        return;
    fi;
    return(sub);
end );

#!
#!
#!  IsIsomorphicSummand(rep1,rep2) only works when rep1 is indecomposable of dimension prime to the field characteristic. 
#!  It returns true if rep1 is isomorphic to a direct summand of rep2, false otherwise. It relies on a result of Benson
#! and Carlson. Written by Peter Webb Feb 2020.
#!
#!


InstallGlobalFunction( IsIsomorphicSummand, function(rep1,rep2)
    local hom, vecs, d;
    if rep1.dimension mod Characteristic(rep1.field) =0 then
        Error("The dimension needs to be prime to the characteristic.");
        return;
    fi;
    if rep1.dimension <> rep2.dimension then
        return false;
    fi;
    hom:=TensorProductRep(DualRep(rep1),rep2);
    vecs:=FixedQuotient(hom);
    d:=Length(vecs);
    vecs:=SafeBaseMat(Concatenation(vecs,FixedPoints@(hom)));
    return Length(vecs) > d;
end );

#!
#! KernelIntersection(rep,rep) . . returns a basis for the intersection of the kernels of all 
#! module homomorphisms A -> B
#! Created by Peter Webb May 8, 2019. Corrected April 18, 2020.
#! This uses a version of SafeNullspaceMat with two arguments, the second being the field.
#!

InstallGlobalFunction( KernelIntersection, function(g,h)
    local transposedmats;
    transposedmats:=List(HomBasis(g,h), TransposedMat);
    if Length(transposedmats)=0 then 
        return(SafeIdentityMat(g.dimension, g.field));
    else
        return SafeNullspaceMat(TransposedMat(Concatenation(transposedmats)),g.field);
    fi;
end );

#
#!PrincipalIdempotent(group,prime) returns a vector in the representation space of 
#!RegularRep(group, GF(prime)) that is the block idempotent of the principal block. Thus if b is 
#!this idempotent and e denotes the vector in the regular representation corresponding to the
#!identity elements, the vector b.e is returned. Spinning it gives a basis for the principal block. 
#!This idempotent is constructed using a formula of B. Kulshammer Arch. Math 56 (1991) 313-319.
#!Code written by Peter Webb February 2020.
#


InstallGlobalFunction( PrincipalIdempotent, function(group,prime)
    local els, primepart, coeffs, pels, qels, qelspos, x, y, position; 
    els:=RightCosetReps(group, Group(Identity(group)));
    primepart:=prime^LogInt(Size(group),prime);
    coeffs:=List(els,x->0);
    pels:=ShallowCopy(els);
    for x in [1..Length(pels)] do
        if RemInt(primepart,Order(pels[x])) <>0 then
            Unbind(pels[x]);
        fi;
    od;
    #pels:=Compacted(pels);
    qels:=ShallowCopy(els);
    qelspos:=List(qels,x->1);
    for x in [1..Length(qels)] do
        if GcdInt(Order(qels[x]),prime)<>1 then
        Unbind(qels[x]);
        qelspos[x]:=0;
        fi;
    od;
    #qels:=Compacted(qels);
    for x in pels do
        for y in qels do
        position:=Position(els, x*y);
        coeffs[position]:=coeffs[position]+qelspos[position];
        od;
    od;
    coeffs:=(Length(qels)*Z(prime)^0)^-1*coeffs;
    return(coeffs);
end );

# -------------------------------- NEW FEATURES IN 1.0.0 ---------------------------------------- #


InstallGlobalFunction( IrreducibleReps, function(rep)
    local L, irreps, gens, i, rho;
    if not IsRepr(rep) then
        Error("for usage, see ?IrreducibleReps");
    fi;
    irreps := IrreducibleRepresentations( rep.group );
    gens := rep.generatorsofgroup;
    L:=[];
    for i in [1..Length( irreps )] do
        rho := GroupHomomorphismByImages(rep.group, Group( List( gens, x -> x^irreps[i] ) ) );
        Add( L, rho );
    od;
    rep.irreps := L;
    return rep;
end );


InstallGlobalFunction( IsRepr, function(rep)
    if IsRecord(rep) then
        if "isRepresentation" in Set(RecNames( rep )) then
            return rep.isRepresentation;
        else
            return false;
        fi;
    else
        return false;
    fi;
end );


InstallGlobalFunction( ReprForG, function(arg)
    # change L by images
    local G, images, rep, arrgim, rho, Ggens;
    arrgim:=arg;
    if not Length(arrgim) in [2..3] or not IsGroup(arrgim[1]) then
        Error("for usage, see ?ReprForG");
    fi;
    G := arrgim[1];
    images := arrgim[2];
    if images = [] then
        Error("Identity group encountered: creating representations is not properly set up for the identity group.");
        return;
    fi;
    if Length(arrgim) = 2 then
        rho := GroupHomomorphismByImages(G, Group( images ) );
    elif Length(arrgim) = 3 then
        Ggens := arrgim[3];
        rho := GroupHomomorphismByImages(G, Group( images ), Ggens, images );
    fi;
    return rho;
end );


InstallGlobalFunction( IrreducibleRepsOfGroup, function(G)
    local L, irreps, gens, i, rho;
    if not IsGroup(G) then
        Error("for usage, see ?IrreducibleRepsOfGroup");
    fi;
    irreps := IrreducibleRepresentations( G );
    gens := GeneratorsOfGroup(G);
    L:=[];
    for i in [1..Length( irreps )] do
        rho := GroupHomomorphismByImages( G, Group( List( gens, x -> x^irreps[i] ) ) );
        Add( L, rho );
    od;
    return L;
end );


InstallGlobalFunction( TensorProductReps, function( rep1, rep2 )

    local gens, mgens, genimages1, genimages2, rho;
    if Source( rep1 ) <> Source( rep1 ) then
        Error("You must have two representations of the same group");
        return;
    fi;
    gens := GeneratorsOfGroup( Source( rep1 ) );
    genimages1 := GeneratorsOfGroup( Images( rep1 ) );
    genimages2 := GeneratorsOfGroup( Images( rep2 ) );
    mgens := List( [1..Length( gens )],
                    i -> TensorProductMatrix( genimages1[i], genimages2[i] ) );
    rho := GroupHomomorphismByImagesNC( Source( rep1 ), Group( mgens ) );
    return rho;
end );

InstallGlobalFunction( DiagonalRep, function( rep )
    local rho;
    rho := REPN_ComputeUsingSerre(rep);
    if IsBound\.( rho, RNamObj( "diagonal_rep" ) ) then
        return rho.diagonal_rep;
    fi;
end );

# -------------------------------- END NEW FEATURES ---------------------------------------- #

GroupRepOps:=rec(Decompose:=DecomposeGroupRep,
                            SubmoduleRep:=SubmoduleGroupRep,
                            QuotientRep:=QuotientGroupRep,
                            Spin:=GroupSpin,
                            CoSpin:=GroupCoSpin,
                            HomBasis:=GroupHomBasis,
                            DecomposeSubmodule:=GroupDecomposeSubmodule,
                            SumOfImages:=GroupSumOfImages,
                            TensorProductRep:=TensorProductGroupRep,
                            DirectSumRep:=DirectSumGroupRep			
                            );