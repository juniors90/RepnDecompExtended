DeclareGlobalFunction( "Rep" );

#! TensorProductMatrix( A, B ) . . . Tensor product of two matrices
#!
#! The GAP command KroneckerProduct seems inexplicably slow and in tests on
#! two 60 x 60 matrices takes about twice as long as the following code.

DeclareGlobalFunction( "TensorProductMatrix" );

DeclareGlobalFunction( "TensorProductReps" );

DeclareGlobalFunction( "IrreducibleRepsOfGroup" );
