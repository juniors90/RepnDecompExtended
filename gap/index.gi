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
    local gens, mgens, gen1, gen2, newrho, rep;
    gens  := rep1.generatorsofgroup;
    gen1  := rep1.genimages;
    gen2  := rep2.genimages;
    mgens := List( [1..Length( gens )],
                    i -> TensorProductMatrix( gen1[i], gen2[i] ) );
    newrho   := GroupHomomorphismByImagesNC( rep1.group , Group( mgens ) );
    rep   := rec(
                group             := rep1.group,
                generatorsofgroup := gens,        # new feature
                rho               := newrho,      # new feature
                irreps            := rep1.irreps, # new feature
                genimages         := mgens,
                isRepresentation  := true,
                isIrreps          := newrho in rho.irreps,
                dimension         := Length(mgens[1]),
                operations        := GroupRepOps,
    );
    return rep;
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


InstallGlobalFunction( DiagonalRep, function( rep )
    local rho;
    rho := REPN_ComputeUsingSerre(rep.rho);
    if IsBound\.( rho, RNamObj( "diagonal_rep" ) ) then
        rep.basis             := rho.basis;
        rep.centralizer_basis := rho.centralizer_basis;
        rep.decomposition     := rho.decomposition;
        rep.diagonal_rep      := rho.diagonal_rep; 
    fi;
    return rep;
end );