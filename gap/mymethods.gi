InstallGlobalFunction( Rep, function(arg)
    # change L by images
    local G, images, rep, arrgh, rho, irreps;
    arrgh:=arg;
    G      := arrgh[1];
    images := arrgh[2];
    if not Length(arrgh) in [2..3] or not IsGroup(G) or IsMatrix(images) then
        Error("for usage, see ?Rep");
    fi;
    if images = [] then
        Error("Identity group encountered: creating representations is not properly set up for the identity group.");
        return;    
    fi;
    if   Length(arrgh) = 2 then
        rho := GroupHomomorphismByImages(G, Group( images ) );
    fi;
    irreps := IrreducibleRepsOfGroup( G );
    rep:=rec(
            group             := G,
            generatorsofgroup := GeneratorsOfGroup( G ),      # new feature
            rho               := rho,                         # new feature
            irreps            := irreps, # new feature
            genimages         := images,
            isRepresentation  := true,
            isIrreps          := rho in irreps,
            dimension         := Length(images[1]),
            operations        := GroupRepOps,
        );
    if Length(arrgh) = 2 then
        rep.field := Field(images[1][1]);
    else
        rep.field := arrgh[1];
    fi;
    return rep;
end );

#! TensorProductMatrix( A, B ) . . . Tensor product of two matrices
#!
#! The GAP command KroneckerProduct seems inexplicably slow and in tests on
#! two 60 x 60 matrices takes about twice as long as the following code.


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


InstallGlobalFunction( TensorProductReps, function( rep1, rep2 )
    local gens, mgens, genimages1, genimages2, rho;
    if Source( rep1 ) <> Source( rep1 ) then
        Error("You must have two representations of the same group");
    fi;
    gens       := GeneratorsOfGroup( Source( rep1 ) );
    genimages1 := GeneratorsOfGroup( Images( rep1 ) );
    genimages2 := GeneratorsOfGroup( Images( rep2 ) );
    mgens      := List( [1..Length( gens )],
                   i -> TensorProductMatrix( genimages1[i], genimages2[i] ) );
    rho        := GroupHomomorphismByImagesNC( Source( rep1 ), Group( mgens ) );
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


#InstallGlobalFunction( DiagonalRep, function( rep )
#    local rho;
#    rho := REPN_ComputeUsingSerre(rep);
#    if IsBound\.( rho, RNamObj( "diagonal_rep" ) ) then
#        return rho.diagonal_rep;
#    fi;
#    return fail;
#end );