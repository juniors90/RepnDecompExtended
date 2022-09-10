
#  ***************************************************************************
#                       Catreps
#           Copyright (C) 2008 Peter Webb
#       Copyright (C) 2011 Peter Webb, Fan Zhang
#       Copyright (C) 2020 Moriah Elkin
#       Copyright (C) 2022 Ferreira Juan David
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#  The overall structure of the catreps package was designed and most if it
#  written by Peter Webb <webb@math.umn.edu>, who is also the maintainer. 
#  Contributions were made by Dan Christensen,
#  Fan Zhang and Moriah Elkin.
#  ***************************************************************************

#!
#! SupportOfMorphism(m) returns the list of positions at which the 
#! list m is defined.
#!

DeclareGlobalFunction( "SupportOfMorphism" );

#!
#! IdentityMorphism(l) returns the identity morphism on the set l.
#! A morphism is stored as a list of the images of its values on a set
#! of numbers, which form its domain of definition, and are taken to be
#! an object in a concrete category. At elements of other objects the
#! morphism will be undefined.
#!

DeclareGlobalFunction( "IdentityMorphism" );

#!
#! Composition(f,g) returns the composition of the functions f and g, expressed
#! as lists of their values.
#!

DeclareGlobalFunction( "Composition" );

#!
#! IsComposable(f,g) returns true if the functions f and g, expressed
#! as lists of their values, can be composed and false otherwise.
#!

DeclareGlobalFunction( "IsComposable" );

#!
#! Objects(cat) returns the objects of the concrete category cat, as a list
#! of sets. At the moment it will not work unless for every object there
#! is at least one generator morphism whose support is that object.
#!

DeclareGlobalFunction( "Objects" );

#!
#! Origin(cat,m) returns the position in cat.objects of the domain of the
#! morphism m.
#!

DeclareGlobalFunction( "Origin" );


#!
#! Terminus(cat,m) returns the position in cat.objects of the domain of the
#! morphism m.
#!

DeclareGlobalFunction( "Terminus" );



DeclareGlobalFunction( "Domains" );

#!
#! ConcreteCategory(list of functions, (list of sets))
#! There are optionally one or two arguments. The first is a list of generating
#! functions of the category, the second is a list of the objects of the category.
#! The function starts a record for a concrete category.
#! If there is only one argument, the objects are taken to be the domains of the
#! generator morphisms, so for every object there should be
#! at least one generator morphism whose support is that object.
#! It could be the identity morphisms, but doesn't have to be.
#!
#! Written by Peter Webb 2008, Moriah Elkin 2018.
#!

DeclareGlobalFunction( "ConcreteCategory" );

#!
#! EmptyMat(r,s) returns an r x s matrix in which each entry is [].
#!

DeclareGlobalFunction( "EmptyMat" );


#!
#! Morphisms(cat) returns an l x l matrix, where l is the number of objects
#! in the category cat, and where the i,j entry is a list of the
#! morphisms from object i to
#! object j.
#!

DeclareGlobalFunction( "Morphisms" );

#!
#! CatRep(category, matrices, field)
#! This function creates a representation of a category where the generators 
#! are represented by the matrices in the list.
#!

DeclareGlobalFunction( "CatRep" );

#!
#! YonedaRep(category, object number, field)
#! This returns the representation of the category on the subspace of the
#! category algebra spanned by the morphisms whose domain is the specified
#! object. The representation is projective, by Yoneda's lemma.
#!


DeclareGlobalFunction( "YonedaRep" );

#!
#! YonedaDualRep(category, object number, field)
#! This returns the representation of the category on the dual of the subspace of the
#! category algebra spanned by the morphisms whose codomain is the specified
#! object. The representation is injective, because its dual is projective, by Yoneda's lemma.
#!


DeclareGlobalFunction( "YonedaDualRep" );

#!
#! ZeroCatRep(cat,field) returns the zero representation of the category.
#!

DeclareGlobalFunction( "ZeroCatRep" );

#!
#!
#! ConstantRep(cat,field) returns the constant (or trivial) representation of the category.
#!
#!

DeclareGlobalFunction( "ConstantRep" );

#!
#! TensorProductCatRep(rep,rep) . . . Kronecker product of two representations
#! Called by TensorProductRep(rep,rep) if rep is a category representation
#!

DeclareGlobalFunction( "TensorProductCatRep" );

#!
#! ExtractHom(vec, dimension vector 1, dimension vector 2)
#! takes a vector and returns it repackaged as a list of dimvec1[i] by dimvec2[i] matrices.
#!

DeclareGlobalFunction( "ExtractHom" );

#!
#! HomBasis(M, N) returns a basis for the space of 
#! kC-module homomorphisms M -> N. Each element of the output list is a
#! list of matrices, which the ith matrix is a map M[i] -> N[i].
#!
#! The algorithm used finds homomorphisms X so that Xg-gX=0 for all category generators
#! g. The code is an elaboration by Peter Webb (October 2008) of code written 
#! for group representations
#! by Dan Christensen in August 2007.
#!
#! CatHomBasis is called by HomBasis(M, N) when M is a
#! category representation.
#!


DeclareGlobalFunction( "CatHomBasis" );

#!
#! SumOfImages(rep,rep) . . returns a list of bases which give the sum of images of all 
#! module homomorphisms A -> B. Term i in the list is a basis for the subspace of the 
#! value of representation B at object i which is part of the sum of images.
#!
#! CatSumOfImages is called by SumOfImages(rep,rep) when rep is a
#! category representation.
#!

DeclareGlobalFunction( "CatSumOfImages" );

#!
#!
#! SubmoduleRep(rep, list of lists of vecs) . .  returns the representation which gives
#! the action on the submodule spanned at each object by the corresponding
#! list of vectors. Each list of vectors must be a basis. This is not checked.
#!
#! SubmoduleCatRep is called by SubmoduleRep(rep, list of lists of vecs) when rep is a
#! category representation.
#!
#!

DeclareGlobalFunction( "SubmoduleCatRep" );

#!
#! QuotientRep(rep, basis structure) . . . returns the representation giving the
#! action on the quotient by an invariant subspace.
#!
#! At the moment this does not work if the basis structure is for the zero subspace.
#!
#! QuotientCatRep is called by QuotientRep(rep, basis structure) when rep is a
#! category representation.
#!

DeclareGlobalFunction( "QuotientCatRep" );

#!
#! DecomposeSubmodule(repn, basis structure) . . probably returns a list of two basis 
#! structures for summands of the submodule of repn spanned by 
#! the given basis structure, if the module is decomposable, and
#! if the module is indecomposable it returns a list whose only element is the
#! given basis.
#!
#! CatDecomposeSubmodule is called by DecomposeSubmodule(repn, basis structure) when repn is a
#! category representation.
#!
#! This code was developed by Peter Webb from code for group representations
#! written by Bryan Simpkins and Robert Hank, University of
#! Minnesota, April 2006, related to an approach taken by Michael Smith in
#! code written in 1993, and tidied up by Dan Christensen, University of Western Ontario, Aug 2007.
#!

DeclareGlobalFunction( "CatDecomposeSubmodule" );

#!
#! Decompose(rep) . . returns a list of basis structures of direct summands of rep which
#! are probably indecomposable.
#!
#! DecomposeCatRep is called by Decompose(rep) when rep is a category representation.
#!
#! This code was developed by Peter Webb from code for group representations
#! written by Bryan Simpkins and Robert Hank, University of
#! Minnesota, April 2006, related to an approach taken by Michael Smith in
#! code written in 1993, and tidied up by Dan Christensen, University of Western Ontario, Aug 2007.
#!

DeclareGlobalFunction( "DecomposeCatRep" );

#!
#! Spin(rep, list of lists of vectors)
#! returns a basis for the subrepresentation generated by the vectors in the lists.
#! There is one entry in the list (of lists of vectors) for each object in the category,
#! and it is a list of vectors which all belong to the representation space at that object.
#!
#! CatSpin is called by Spin(rep, list of lists of vectors) when rep is a
#! category representation.
#!
#! This is code which was developed from code written by
#! Fan Zhang (University of Minnesota), March 2011.
#!

DeclareGlobalFunction( "CatSpin" );

#!
#! OrthogonalComplement(list of vectors, dim, field)
#! returns a basis for the orthgonal complement of the space spanned by the
#! list of vectors, in a space of dimension dim (put there in case the list
#! of vectors is empty).
#!

DeclareGlobalFunction( "OrthogonalComplement" );

#!
#! CoSpin(rep, list of lists of vectors)
#! returns a basis for the largest subrepresentation contained in the
#! subspaces spanned by the vectors in the lists.
#! There is one entry in the list (of lists of vectors) for each object in the category,
#! and it is a list of vectors which all belong to the representation space at that object.
#!
#! CatCoSpin is called by CoSpin(rep, list of lists of vectors) when rep is a
#! category representation.
#!

DeclareGlobalFunction( "CatCoSpin" );

#!
#! EndomorphismGroup(cat, obj) returns the group of the endomorphisms of obj in
#! cat, in permutation form. It assumes every endomorphism is invertible and that the generators of the 
#! endomorphism group appear among the generators of the category.
#!
#! Written by Moriah Elkin July 2018.
#!

DeclareGlobalFunction( "EndomorphismGroup" );

#!
#! EndomorphismGroups(cat) returns a list containing the endomorphism group for
#! each object in the category cat.
#!
#! Written by Moriah Elkin July 2018.
#!

DeclareGlobalFunction( "EndomorphismGroups" );

#!
#! FI(n) and FI2(n) interchangeably return a record for the category FI with
#! objects 0...n. O is represented by the first object ( [1] ), and its morphisms
#! correspond to the first element in every object, which is otherwise ignored
#! (first elements only map to first elements). The category FI is the category of finite sets with
#! injective maps that has featured in the theory of representation stability.
#!
#! Written by Moriah Elkin July 2018.
#!

DeclareGlobalFunction( "FI" );
DeclareGlobalFunction( "FI2" );

#!
#! DirectSumCatRep(rep1, rep2) returns the representation that is the direct
#! sum of the representations rep1 and rep2 of the category rep1.category.
#!
#! Written by Moriah Elkin July 2018.
#!
#! DirectSumCatRep is called by DirectSumRep(rep1,rep2) when rep1 and rep2 are
#! group representations.
#!

DeclareGlobalFunction( "DirectSumCatRep" );

#!
#! GeneratorDomains(cat) returns a l x l matrix, where l is the number of objects
#! in the category cat, and where the i,j entry is a list of the indices of morphisms
#! from object i to object j.
#!
#! Written by Moriah Elkin August 2018.
#!

DeclareGlobalFunction( "GeneratorDomains" );

#!
#! MorphismsRep(rep) returns a l x l matrix, where l is the number of objects
#! in the category of rep, and where the i,j entry is a list of the matrices of
#! morphisms from object i to object j.
#!
#! Written by Moriah Elkin (August 2018), based on code for categories written
#! by Peter Webb.
#!

DeclareGlobalFunction( "MorphismsRep" );

#!
#! SubConstant(rep) returns a list of bases for the largest subconstant
#! representation of rep. The algorithm finds the common kernel of all differences
#! of morphisms from a domain to all codomains.
#!
#! Written by Moriah Elkin August 2018.
#!

DeclareGlobalFunction( "SubConstant" );

#!
#! Evaluation(rep, obj) returns the representation of the endomorphism
#! group of object obj on the value of the representation rep at obj. (In FI, obj is
#! not the mathematical object in FI, but rather the object number in the category).
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "Evaluation" );

#!
#! SpinFixedPts(rep, obj) returns a list of lists of bases (lists of vectors),
#! where list of bases i is the spin of the fixed points of the group representation
#! that is the i'th summand of the evaluation of rep at object obj; item i in the list
#! is an empty list when there are no generators or no fixed points in that summand of
#! the evaluation.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "SpinFixedPts" );

#!
#! SpinBasisVec(rep, obj) goes through each summand in the evaluation of rep at obj,
#! takes the first vector in its basis, and returns a basis for the subfunctor generated
#! by this vector (an empty list if the basis is empty). It returns a list of lists of
#! bases, where the i'th list of bases corresponds to the i'th summand of the evaluation.
#! It does not check that the summand does not have fixed points.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "SpinBasisVec" );

#!
#! IsDirectSum(summands,sum) takes in a list of bases (lists of vectors), summands, and a
#! basis (list of vectors), sum, and returns true if the direct sum of summands is sum. 
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "IsDirectSum" );

#!
#! TestOneFixedPt(rep) tests whether there is at most one fixed point in the summands
#! of the evaluations of rep at each object. It returns true if there is, false if
#! there is not.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "TestOneFixedPt" );

#!
#! DisplayButterflyDims(sn,subGens) takes in the symmetric group of a certain dimension
#! and the generators for a subgroup corresponding to a partition. It displays the
#! ButterflyFactors matrix (without decomposing), and returns the ButterflyFactors.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "DisplayButterflyDims" );

#!
#! ButterflyDimsRep(rep) takes in a representation, displays the ButterflyFactors matrix
#! (without decomposing), and returns the ButterflyFactors.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "ButterflyDimsRep" );

#!
#! npButterflyDimsRep(rep) takes in a representation and returns the ButterflyFactors
#! (without printing anything).
#!
#! Written by Moriah Elkin Spring 2020.
#!

DeclareGlobalFunction( "npButterflyDimsRep" );

#!
#! SafeFixedPoints(rep) finds the fixed points of a representation rep; if rep
#! has dimension 0, it returns an empty list.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "SafeFixedPoints" );

#!
#! SafeDimHom(rep) returns the dimension of the space of module homomorphisms
#! rep1 -> rep2. If rep has dimension 0, it returns 0.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "SafeDimHom" );

#!
#! ExamineButterflyFactors(b,dRec) takes in a ButterflyFactors matrix and a record of
#! possible factors. It displays the original dimension matrix, and then an ordered list of
#! matrices, where the dimensions in each matrix correspond to the dimensions of the
#! homomorphisms between that element in the ButterflyFactors matrix and the factor in the
#! list, or the number of fixed points in that element in the ButterflyFactors matrix.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "ExamineButterflyFactors" );

#!
#! RemoveFromBottom(rep,d) takes representations rep and d.
#! It returns a representation with all images of d removed from the bottom
#! of rep.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "RemoveFromBottom" );



#!
#! FISummandEvalReps(n,obj) returns a list of the representations of the summands of the
#! evaluations of the summands of the Yoneda representation over GF(2) of FI(n) at the
#! mathematical object obj. There will be 3 levels of lists: the overall list will contain
#! a list for each submodule of the Yoneda representation, each of which will contain a list
#! of submodule representations for each object at which it has been evaluated.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "FISummandEvalReps" );

#!
#! ProjSummands(n,obj) returns whether the summands in FISummandEvalReps(n,obj,GF(2))
#! are projective: each column of the output is the evaluation at a different mathematical
#! object, each row is a different summand of the YonedaRep of FI(n) at obj, and each
#! element in the matrix is a list with whether the summands of the evaluation are
#! projective, in order. Representations with no generators display "fail."
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "ProjSummands" );

#!
#! DimSummands(n,obj) returns the dimensions of the summands in
#! FISummandEvalReps(n,obj,GF(2)): each column of the output is the evaluation at a
#! different mathematical object, each row is a different summand of the YonedaRep of FI(n)
#! at obj, and each element in the matrix is a list with the dimensions of the summands of 
#! the evaluation in order.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "DimSummands" );

#!
#! FixedPtDims(n) returns a list of lists of lists of dimensions. The outer list contains a
#! list for SpinFixedPts of the Yoneda Rep of FI(n) at mathematical object 2, at each object
#! in FI(n); and each of those lists contains the dimensions of the corresponding list of
#! bases at each object (an empty list if that summand has no fixed points).
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "FixedPtDims" );

#!
#! FixedPtBases(n) returns a list of lists of bases. The outer list contains a list for
#! SpinFixedPts of the Yoneda Rep of FI(n) at mathematical object 2, at each object in
#! FI(n); and each of those lists contains the corresponding list of bases at each object
#! (taken at the summand that has fixed points).
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "FixedPtBases" );

#!
#! BasisVecDims(n) returns a list of lists of lists of dimensions. The outer list contains a
#! list for SpinBasisVec of the Yoneda Rep of FI(n) at mathematical object 2, at each object
#! in FI(n); and each of those lists contains the dimensions of the corresponding list of
#! bases at each object (where dimension list i contains the dimensions of
#! SpinBasisVec(rep,i)) (an empty list if that summand's basis is empty). I.e.,
#! BasisVecDims[i][j] contains the lengths of the bases generated by spinning the first
#! vector in the basis for the jth summand of the evaluation of the Yoneda rep at object i.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "BasisVecDims" );

#!
#! BasisVecBases(n) returns a list of lists of lists of bases. The outer list contains a 
#! list for SpinBasisVec of the Yoneda Rep of FI(n) at mathematical object 2, at each object
#! in FI(n); and each of those lists contains the corresponding list of bases at each object
#! (where list i contains SpinBasisVec(rep,i)) (an empty list if that summand's basis is
#! empty). I.e., BasisVecBases[i][j] contains the list of bases generated by spinning the
#! first vector in the basis for the jth summand of the evaluation of the Yoneda rep at
#! object i.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "BasisVecBases" );

#!
#! TestYRSummand(basis,eval,summandi) takes in a basis, an evaluation of a Yoneda Rep over
#! GF(2), and the index of the summand of that evaluation that the basis is thought to be
#! equivalent to; it returns whether it is in fact equivalent.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "TestYRSummand" );

#!
#! AllSpinDims(n, evalObj,outputRec) spins all vectors in the evaluation at math object
#! evalObj of the Yoneda representation of FI(n) at math object 2 with respect to that
#! Yoneda representation. If outputRec is false, returns a sorted list of the dimensions of
#! the bases; if it is true, returns a record with each dimension list corresponding to the 
#! list of vectors that generated it.
#!
#! Written by Moriah Elkin Spring 2019.
#!

DeclareGlobalFunction( "AllSpinDims" );

CatRepOps:=rec(Decompose:=DecomposeCatRep,
                            SubmoduleRep:=SubmoduleCatRep,
                            QuotientRep:=QuotientCatRep,
                            Spin:=CatSpin,
                            CoSpin:=CatCoSpin,
                            HomBasis:=CatHomBasis,
                            DecomposeSubmodule:=CatDecomposeSubmodule,
                            SumOfImages:=CatSumOfImages,
                            TensorProductRep:=TensorProductCatRep,
                            DirectSumRep:=DirectSumCatRep
                            );