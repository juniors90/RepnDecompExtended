InstallGlobalFunction( PermutationRep, function( group, perms, field )
    local d;
    d := Maximum( List( perms, LargestMovedPointPerm ) );
    return Rep( group, List( perms, x -> PermutationMatrix( x, d, field ) ) );
end );

InstallGlobalFunction( PermutationRepOnCosets, function(G, H, F)
    local n, L, M, phi, g, i, j, R, glg, k;
    L := [];
    R := [];
    L := RightCosets( G, H );
    n := Length(L);
    for k in [1..Length( GeneratorsOfGroup( G ) )] do
        g := GeneratorsOfGroup( G )[k];
        M := NullMat( n, n, F);
        for i in [1..n] do
            j := 0;
            repeat j := j + 1;
            until L[i]*g = L[j];
            M[i][j] := Z( Size( F ) )^0;
        od;
        R[k] := M;
    od;
    phi := Rep(G,R,F);
    return phi;
end );



#! SubmoduleGroupRep is called by SubmoduleRep(rep, list of vecs) when rep is a
#! group representation.

InstallGlobalFunction( SubmoduleRep, function( rep, v )
    return rep.operations.SubmoduleRep(rep, v);
end );


InstallGlobalFunction( RemoveFromTop, function ( rep1, rep2 )
    local  newrep;
    newrep := SubmoduleRep( rep1, KernelIntersection( rep1, rep2 ) );
    return newrep;
end );