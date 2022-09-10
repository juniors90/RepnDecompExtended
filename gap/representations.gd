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

#! @Chapter Commands that return a representation
#!
#! Commands that return a representation.
#!
#! @Section Efficient summing over groups

#! @Arguments group, matrices, field(optional)

#! @Returns A representation in which the generators of the group
#! act by the matrices in the list. This and the following commands
#! are the basic way to input a representation. The function takes
#! two or three arguments. If only two arguments are given the field
#! is taken to be the field of the leading entry of the first matrix,
#! which might be incorrect.

#! @Description This function creates a representation provided
#! the group is not the identity group and the representation is
#! not the zero representation.
#!
#! Modification by Craig Corsi and Peter Webb (University of Minnesota)
#! 2016 to allow the optional field argument.
DeclareGlobalFunction( "Rep" );

#! @Arguments group, perms, field

#! @Returns The representation on the permutation module determined
#! by the permutation representation in which the group generators
#! act as the permutations in the `perms` list.

#! @Description `PermutationRep(group, permns, field)` returns
#! a permutation representation in which the group generators act via
#! the permutations in the list.
DeclareGlobalFunction( "PermutationRep" );

#! @Arguments G, H, F

#! @Returns The representation which is the permutation module
#! determined by the right cosets of the subgroup `H` in 
#! the group `G` over the field `F`.
#!
#! Modification by Craig Corsi (University of Minnesota)
#! 2016 to make sure the field is correctly recorded.
DeclareGlobalFunction( "PermutationRepOnCosets" );

#! PermGroupToRep(group, field). . returns the representation on the defining permutation representation.

#! RegularRep(group, field) returns the group ring as a representation.
#! FreeRep(group, field, n) returns the direct sum of n copies of the group ring as a representation.
#! TrivialRep(group, field). . returns a trivial representation.
#! TrivialRepN(group, field, n).. returns a trivial action on an n-dimensional vector space.
#! ZeroGroupRep(group, field). . returns the zero representation.
#! DualRep(rep). . returns a representation which is the dual of rep.
#! DirectSumRep(rep, rep). . returns a representation which is the direct sum of the two representations.
#! TensorProductRep(rep, rep). . returns a representation which is the tensor product of the two representations.

#! SubmoduleRep(rep, list of vecs). . returns the representation on the submodule spanned by the list of vectors. The list of vectors must be a basis for an invariant subspace - this is not checked.


#! QuotientRep(rep, list of vecs). . returns the representation on the quotient module by the submodule spanned by the list of vectors. The list of vectors must be a basis for an invariant subspace - this is not checked.
#! SectionRep(rep, list of vectors, list of vectors). .returns the representation on the quotient of the submodule spanned by the first list of vectors, by the submodule spanned by the second list of vectors. The two lists should be independent sets, and should be bases for submodules of rep, the second submodule contained in the first. This is not checked.
#! RestrictedRep(group, subgroup, rep of group). . returns the restriction of the representation.
#! InducedRep(group, subgroup, rep of subgroup). . returns the induced representation.
#! SymmetricPowerRep(rep, n). . returns the representation on the nth symmetric power of the given representation. This code was written by Brad Froehle.
#! CompositionFactorsRep(rep). . returns a list of representations which are the composition factors of rep. This function calls the meataxe routine MTX.CompositionFactors which is already implemented in GAP.
#! ProjectiveFreeCore(rep). . returns a representation isomorphic to a maximal summand of rep that has no projective summand. It is only guaranteed to work when the group is a p-group in characteristic p. In other cases it may give a correct answer, and if it does not then an error message is returned. The algorithm uses fewer resources than ProjectiveFreeSummand.
#! ProjectiveFreeSummand(rep). . returns a representation isomorphic to a maximal summand of rep that has no projective summand. The algorithm decomposes rep and tests summands for projectivity. This is computationally expensive.
#! FixedPointRep(rep, subgroup of rep.group). .returns the representation of the normalizer of the subgroup on the fixed points of the subgroup.
#! BrauerRep(rep, p-subgroup of rep.group). .returns the representation of the normalizer of the p-subgroup on the Brauer quotient at that subgroup (quotient of the fixed points by the images of trace maps from proper subgroups).
#! RemoveFromBottom(rep1, rep2). .returns a representation with the sum of all images of rep2 in rep1 factored out from rep1.

#! @Arguments rho, tau
#!
#! @Returns a representation that is the common kernel
#! of all homomorphisms from `rho` to `tau`.
#!
#! @Description Written by Peter Webb Spring 2020 
DeclareGlobalFunction( "RemoveFromTop" );