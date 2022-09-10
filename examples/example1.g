# ***************************************************************************
#                                   Reps
#       Copyright (C) 2005 Peter Webb 
#       Copyright (C) 2006 Peter Webb, Robert Hank, Bryan Simpkins 
#       Copyright (C) 2007 Peter Webb, Dan Christensen, Brad Froehle
#       Copyright (C) 2020 Peter Webb, Moriah Elkin
#       Copyright (C) 2022 Peter Webb, Ferreira Juan David
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#  The overall structure of the reps package was designed and most of it
#  written by Peter Webb <webb@math.umn.edu>, who is also the maintainer. 
#  Contributions were made by
#  Dan Christensen, Roland Loetscher, Robert Hank, Bryan Simpkins,
#  Brad Froehle and others.
# ***************************************************************************
#
# LoadPackage("RepnDecompExtended", "0", false);
#! @BeginChunk Example_CanonicalDecomposition
#! @BeginExample
#! This is the Symmetric Group S_3
G := SymmetricGroup( 3 );
#! Sym( [ 1 .. 3 ] )
m := [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ];
rep := Rep(G,m);
#! rec( dimension := 1,
#!        field := Rationals,
#!        generatorsofgroup := [ (1,2,3), (1,2) ],
#!        genimages := [ [ [ 1 ] ], [ [ 1 ] ] ],
#!        group := Sym( [ 1 .. 3 ] ),
#!  irreps := [
#!    [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ],
#!    [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ],
#!    [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ],
#!                          [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] ],
#!  isIrreps := true,
#!  isRepresentation := true,
#!  operations := rec(  ),
#!  rho := [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ] )
# Calculamos las representaciones irrreducibles para el
# grupo simétrico $S_3$.
# Primera parte: calculamos la Suma directa de las
# representationes irreducibles
triv_sum_sgn_sum_std := DirectSumOfRepresentations(rep.irreps);
#! [ (1,2,3), (1,2) ] -> [ [ [ 1, 0, 0, 0 ],
#!                           [ 0, 1, 0, 0 ],
#!                           [ 0, 0, E(3), 0 ],
#!                           [ 0, 0, 0, E(3)^2 ] ],
#!                         [ [ 1, 0, 0, 0 ],
#!                           [ 0, -1, 0, 0 ],
#!                           [ 0, 0, 0, E(3)^2 ],
#!                           [ 0, 0, E(3), 0 ] ] ]
IrreducibleDecomposition(triv_sum_sgn_sum_std);
#! [ rec( basis := [ [ 1, 0, 0, 0 ] ] ),
#!   rec( basis := [ [ 0, 1, 0, 0 ] ] ),
#!   rec( basis := [ [ 0, 0, 1, 0 ],
#!                   [ 0, 0, 0, 1 ] ] ) ]
# Calculamos una la representación de $S_3$. dada por el
# producto tensorial de la representación estandar por
# la representación estandar.
std_tensor_std := TensorProductReps(rep, rep);
#! [ (1,2,3), (1,2) ] -> [ [ [ E(3)^2, 0, 0, 0 ],
#!                           [ 0, 1, 0, 0 ],
#!                           [ 0, 0, 1, 0 ],
#!                           [ 0, 0, 0, E(3) ] ],
#!                         [ [ 0, 0, 0, E(3) ],
#!                           [ 0, 0, 1, 0 ],
#!                           [ 0, 1, 0, 0 ],
#!                           [ E(3)^2, 0, 0, 0 ] ] ]
# Calculamos base, representación diagonal, decomposition
# y base centralizadora de la representación de S3 dada
# por el producto tensorial de la representación estandar
# por la representación estandar.
IrreducibleDecomposition(std_tensor_std.rho);
#! [ rec( basis := [ [ 0, 1, 1, 0 ] ] ),
#!   rec( basis := [ [ 0, 1, -1, 0 ] ] ),
#!   rec( basis := [ [ 0, 0, 0, 1 ],
#!                   [ 1, 0, 0, 0 ] ] ) ]
std_tensor_std_diag_rep :=  DiagonalRep(std_tensor_std);
#! [ (1,2,3), (1,2) ] -> [ [ [ 1, 0, 0, 0 ],
#!                           [ 0, 1, 0, 0 ],
#!                           [ 0, 0, E(3), 0 ],
#!                           [ 0, 0, 0, E(3)^2 ] ],
#!                         [ [ 1, 0, 0, 0 ],
#!                           [ 0, -1, 0, 0 ],
#!                           [ 0, 0, 0, E(3)^2 ],
#!                           [ 0, 0, E(3), 0 ] ] ]
#! [ rec( basis := [ [ 1, 0, 0, 0 ] ] ),
#!   rec( basis := [ [ 0, 1, 0, 0 ] ] ),
#!   rec( basis := [ [ 0, 0, 1, 0 ],
#!                   [ 0, 0, 0, 1 ] ] ) ]
std_tensor_std_diag_rep.diagonal_rep = triv_sum_sgn_sum_std;
#! true
#! @EndExample
#! @EndChunk