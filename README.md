# The GAP package 

[![Build status](https://github.com/juniors90/RepnDecompExtended/actions/workflows/CI.yml/badge.svg)](https://github.com/juniors90/RepnDecompExtended/actions)
[![codecov](https://codecov.io/github/juniors90/RepnDecompExtended/branch/main/graph/badge.svg?token=anQ5vNi9MQ)](https://codecov.io/github/juniors90/RepnDecompExtended)
[![GitHub issues](https://img.shields.io/github/issues/juniors90/RepnDecompExtended)](https://github.com/juniors90/RepnDecompExtended/issues)
[![GitHub forks](https://img.shields.io/github/forks/juniors90/RepnDecompExtended)](https://github.com/juniors90/RepnDecompExtended/network)
[![GitHub stars](https://img.shields.io/github/stars/juniors90/RepnDecompExtended)](https://github.com/juniors90/RepnDecompExtended/stargazers)
[![GitHub contributors](https://img.shields.io/github/contributors/juniors90/RepnDecompExtended?color=green)](https://github.com/juniors90/RepnDecompExtended/graphs/contributors)


```gap
gap> LoadPackage("RepnDecompExtended", "0", false);
true
gap> G := SymmetricGroup( 3 );
Sym( [ 1 .. 3 ] )
gap> m1:=[ [ E(3), 0 ], [ 0, E(3)^2 ] ];;
gap> m2:=[ [ 0,    E(3)^2 ], [ E(3), 0] ];;
gap> m := [ m1, m2 ];;
gap> rep := Rep(G,m);
rec( dimension := 2, field := CF(3), generatorsofgroup := [ (1,2,3), (1,2) ], 
  genimages := [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ], 
  group := Sym( [ 1 .. 3 ] ), 
  irreps := [ [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ], 
      [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ], 
      [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] 
     ], isIrreps := true, isRepresentation := true, operations := rec(  ), 
  rho := [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], 
      [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] )
gap> rep.rho;
[ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ]
gap> rep.irreps;
[ [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ], [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ], 
  [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] ]
gap> rep.irreps[1];
[ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ]
gap> rep.irreps[2];
[ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ]
gap> rep.irreps[3];
[ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ]
gap> triv_sum_sgn_sum_std := DirectSumOfRepresentations(rep.irreps);;
gap> triv_sum_sgn_sum_std;
[ (1,2,3), (1,2) ] -> [ [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, E(3), 0 ], [ 0, 0, 0, E(3)^2 ] ], 
  [ [ 1, 0, 0, 0 ], [ 0, -1, 0, 0 ], [ 0, 0, 0, E(3)^2 ], [ 0, 0, E(3), 0 ] ] ]
gap> Print( IrreducibleDecomposition(triv_sum_sgn_sum_std) );
[ rec(
      basis := [ [ 1, 0, 0, 0 ] ] ), rec(
      basis := [ [ 0, 1, 0, 0 ] ] ), rec(
      basis := [ [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] ] ) ]
gap> std_tensor_std := TensorProductReps(rep, rep);;
gap> std_tensor_std;
rec( dimension := 4, generatorsofgroup := [ (1,2,3), (1,2) ], 
  genimages := [ [ [ E(3)^2, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, E(3) ] ], 
      [ [ 0, 0, 0, E(3) ], [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ E(3)^2, 0, 0, 0 ] ] ], group := Sym( [ 1 .. 3 ] ), 
  irreps := [ [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ], [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ], 
      [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] ], isIrreps := false, 
  isRepresentation := true, operations := rec(  ), 
  rho := [ (1,2,3), (1,2) ] -> [ [ [ E(3)^2, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, E(3) ] ], 
      [ [ 0, 0, 0, E(3) ], [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ E(3)^2, 0, 0, 0 ] ] ] )
gap>  Print( IrreducibleDecomposition(std_tensor_std.rho) );
[ rec(
      basis := [ [ 0, 1, 1, 0 ] ] ), rec(
      basis := [ [ 0, 1, -1, 0 ] ] ), rec(
      basis := [ [ 0, 0, 0, 1 ], [ 1, 0, 0, 0 ] ] ) ]
gap> std_tensor_std_diag_rep :=  DiagonalRep(std_tensor_std);;
gap> std_tensor_std_diag_rep;
rec( basis := [ [ 0, 1, 1, 0 ], [ 0, 1, -1, 0 ], [ 0, 0, 0, 1 ], [ 1, 0, 0, 0 ] ], 
  centralizer_basis := [ [ [ [ 1 ] ], [ [ 0 ] ], [ [ 0, 0 ], [ 0, 0 ] ] ], 
      [ [ [ 0 ] ], [ [ 1 ] ], [ [ 0, 0 ], [ 0, 0 ] ] ], [ [ [ 0 ] ], [ [ 0 ] ], [ [ 1, 0 ], [ 0, 1 ] ] ] ], 
  decomposition := [ [ rec( basis := [ [ 0, 1, 1, 0 ] ] ) ], [ rec( basis := [ [ 0, 1, -1, 0 ] ] ) ], 
      [ rec( basis := [ [ 0, 0, 0, 1 ], [ 1, 0, 0, 0 ] ] ) ] ], diagonal_rep := [ (1,2,3), (1,2) ] -> 
    [ [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, E(3), 0 ], [ 0, 0, 0, E(3)^2 ] ], 
      [ [ 1, 0, 0, 0 ], [ 0, -1, 0, 0 ], [ 0, 0, 0, E(3)^2 ], [ 0, 0, E(3), 0 ] ] ], dimension := 4, 
  generatorsofgroup := [ (1,2,3), (1,2) ], 
  genimages := [ [ [ E(3)^2, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, E(3) ] ], 
      [ [ 0, 0, 0, E(3) ], [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ E(3)^2, 0, 0, 0 ] ] ], group := Sym( [ 1 .. 3 ] ), 
  irreps := [ [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ 1 ] ] ], [ (1,2,3), (1,2) ] -> [ [ [ 1 ] ], [ [ -1 ] ] ], 
      [ (1,2,3), (1,2) ] -> [ [ [ E(3), 0 ], [ 0, E(3)^2 ] ], [ [ 0, E(3)^2 ], [ E(3), 0 ] ] ] ], isIrreps := false, 
  isRepresentation := true, operations := rec(  ), 
  rho := [ (1,2,3), (1,2) ] -> [ [ [ E(3)^2, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, E(3) ] ], 
      [ [ 0, 0, 0, E(3) ], [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ E(3)^2, 0, 0, 0 ] ] ] )
gap> Print( IrreducibleDecomposition( std_tensor_std_diag_rep.diagonal_rep ) );
[ rec(
      basis := [ [ 1, 0, 0, 0 ] ] ), rec(
      basis := [ [ 0, 1, 0, 0 ] ] ), rec(
      basis := [ [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] ] ) ]
gap> std_tensor_std_diag_rep.diagonal_rep = triv_sum_sgn_sum_std;
true
```

## Contact

TODO: add info on how to contact you and/or how to report issues with your
package

## License

TODO: Provide information on the license of your package. A license is
important as it determines who has a right to distribute your package. The
"default" license to consider is GNU General Public License v2 or later, as
that is the license of GAP itself.
