# The GAP package 

```
LoadPackage("RepnDecompExtended", "0", false); # Cargamos el paquete
G := SymmetricGroup( 3 ); # Definimos el grupo simétrico S3
# Calculamos las representaciones irrreducibles para el
# grupo simétrico
irreps := IrreducibleRepsOfGroup(G);
# Definimos las variables triv3, sgn3 y std3 que contenga a la 
# representación trivial, signo y estandar del grupo simétrico S3
triv3 := irreps[1]; sgn3 := irreps[2]; std3:=irreps[3];
# Calculamos una la representación de S3 dada por el producto
# tensorial de la representación estandar por la representación estandar
std3_tensor_std3 := TensorProductReps(std3, std3);
# Calculamos base, representación diagonal, decomposition y base
# centralizadora de la representación de S3 dada por el producto
# tensorial de la representación estandar por la representación estandar.
irred_decomp := REPN_ComputeUsingSerre(std3_tensor_std3);
# # Guardamo en una variable la representación diagonal
std3_tensor_std3_diagonal_rep := irred_decomp.diagonal_rep;
triv3_sum_sgn3_sum_std3 := DirectSumOfRepresentations([triv3, sgn3, std3]);
irred_decomp_triv3_sum_sgn3_sum_std3 := IrreducibleDecomposition(triv3_sum_sgn3_sum_std3);
std3_tensor_std3_diagonal_rep=triv3_sum_sgn3_sum_std3;
```





LoadPackage("RepnDecompExtended", "0", false);
G := SymmetricGroup( 3 );;
irreps := IrreducibleRepsOfGroup(G);;
triv_sum_sgn_sum_std := DirectSumOfRepresentations(irreps);;
IrreducibleDecomposition(triv_sum_sgn_sum_std);;
triv := irreps[1];;
sgn := irreps[2];;
std := irreps[3];;
std_tensor_std := TensorProductReps(std, std);;
IrreducibleDecomposition(std_tensor_std);;
std_tensor_std_diag_rep := DiagonalRep(std_tensor_std);;
IrreducibleDecomposition(std_tensor_std_diag_rep);;
std_tensor_std_diag_rep = triv_sum_sgn_sum_std;


## Contact

TODO: add info on how to contact you and/or how to report issues with your
package

## License

TODO: Provide information on the license of your package. A license is
important as it determines who has a right to distribute your package. The
"default" license to consider is GNU General Public License v2 or later, as
that is the license of GAP itself.
