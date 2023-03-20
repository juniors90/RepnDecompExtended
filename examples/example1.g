# ***************************************************************************
#                                   Reps
#       Copyright (C) 2022 Ferreira Juan David
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
# ***************************************************************************
#
#! @BeginChunk Example_CanonicalDecomposition
#! @BeginExample
LoadPackage("RepnDecompExtended", "0", false);
#! true
G := SymmetricGroup( 3 );
#! Sym( [ 1 .. 3 ] )
m1:=[ [ E(3), 0 ], [ 0, E(3)^2 ] ];;
m2:=[ [ 0,    E(3)^2 ], [ E(3), 0] ];;
m := [ m1, m2 ];;
rep := Rep(G,m);
#! rec( dimension := 2, field := CF(3), generatorsofgroup := [ (1,2,3), (1,2) ], 
#!   genimages := [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ], 
#!   group := Sym( [ 1 .. 3 ] ), 
#!   irreps := [ [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ], 
#!       [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ], 
#!       [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] 
#!      ], isIrreps := true, isRepresentation := true, operations := rec(  ), 
#!   rho := [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], 
#!       [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] )
rep.rho;
#! [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ]
rep.irreps;
#! [ [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ], [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ], 
#!   [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] ]
rep.irreps[1];
#! [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ]
rep.irreps[2];
#! [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ]
rep.irreps[3];
#! [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ]
triv_sum_sgn_sum_std := DirectSumOfRepresentations(rep.irreps);;
triv_sum_sgn_sum_std;
#! [ (1,2,3), (1,2) ] -> [ [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, E(3), 0 ], [ 0, 0, 0, E(3)^2 ] ], 
#!   [ [ 1, 0, 0, 0 ], [ 0, -1, 0, 0 ], [ 0, 0, 0, E(3)^2 ], [ 0, 0, E(3), 0 ] ] ]
Print( IrreducibleDecomposition(triv_sum_sgn_sum_std) );
#! [ rec(
#!       basis := [ [ 1, 0, 0, 0 ] ] ), rec(
#!       basis := [ [ 0, 1, 0, 0 ] ] ), rec(
#!       basis := [ [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] ] ) ]
std_tensor_std := TensorProductReps(rep, rep);;
std_tensor_std;
#! rec( dimension := 4, generatorsofgroup := [ (1,2,3), (1,2) ], 
#!   genimages := [ [ [ E(3)^2, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, E(3) ] ], 
#!       [ [ 0, 0, 0, E(3) ], [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ E(3)^2, 0, 0, 0 ] ] ], group := Sym( [ 1 .. 3 ] ), 
#!   irreps := [ [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ], [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ], 
#!       [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] ], isIrreps := false, 
#!   isRepresentation := true, operations := rec(  ), 
#!   rho := [ (1,2,3), (1,2) ] -> [ [ [ E(3)^2, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, E(3) ] ], 
#!       [ [ 0, 0, 0, E(3) ], [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ E(3)^2, 0, 0, 0 ] ] ] )
Print( IrreducibleDecomposition(std_tensor_std.rho) );
#! [ rec(
#!       basis := [ [ 0, 1, 1, 0 ] ] ), rec(
#!       basis := [ [ 0, 1, -1, 0 ] ] ), rec(
#!       basis := [ [ 0, 0, 0, 1 ], [ 1, 0, 0, 0 ] ] ) ]
std_tensor_std_diag_rep :=  DiagonalRep(std_tensor_std);;
std_tensor_std_diag_rep;
#! rec( basis := [ [ 0, 1, 1, 0 ], [ 0, 1, -1, 0 ], [ 0, 0, 0, 1 ], [ 1, 0, 0, 0 ] ], 
#!   centralizer_basis := [ [ [ [ 1 ] ], [ [ 0 ] ], [ [ 0, 0 ], [ 0, 0 ] ] ], 
#!       [ [ [ 0 ] ], [ [ 1 ] ], [ [ 0, 0 ], [ 0, 0 ] ] ], [ [ [ 0 ] ], [ [ 0 ] ], [ [ 1, 0 ], [ 0, 1 ] ] ] ], 
#!   decomposition := [ [ rec( basis := [ [ 0, 1, 1, 0 ] ] ) ], [ rec( basis := [ [ 0, 1, -1, 0 ] ] ) ], 
#!       [ rec( basis := [ [ 0, 0, 0, 1 ], [ 1, 0, 0, 0 ] ] ) ] ], diagonal_rep := [ (1,2,3), (1,2) ] -> 
#!     [ [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, E(3), 0 ], [ 0, 0, 0, E(3)^2 ] ], 
#!       [ [ 1, 0, 0, 0 ], [ 0, -1, 0, 0 ], [ 0, 0, 0, E(3)^2 ], [ 0, 0, E(3), 0 ] ] ], dimension := 4, 
#!   generatorsofgroup := [ (1,2,3), (1,2) ], 
#!   genimages := [ [ [ E(3)^2, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, E(3) ] ], 
#!       [ [ 0, 0, 0, E(3) ], [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ E(3)^2, 0, 0, 0 ] ] ], group := Sym( [ 1 .. 3 ] ), 
#!   irreps := [ [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ], [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ], 
#!       [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] ], isIrreps := false, 
#!   isRepresentation := true, operations := rec(  ), 
#!   rho := [ (1,2,3), (1,2) ] -> [ [ [ E(3)^2, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, E(3) ] ], 
#!       [ [ 0, 0, 0, E(3) ], [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ E(3)^2, 0, 0, 0 ] ] ] )
Print( IrreducibleDecomposition( std_tensor_std_diag_rep.diagonal_rep ) );
#! [ rec(
#!       basis := [ [ 1, 0, 0, 0 ] ] ), rec(
#!       basis := [ [ 0, 1, 0, 0 ] ] ), rec(
#!       basis := [ [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] ] ) ]
std_tensor_std_diag_rep.diagonal_rep = triv_sum_sgn_sum_std;
#! true
#! @EndExample
#! @EndChunk
