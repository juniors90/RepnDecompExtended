
FixedPoints@ := function(rep)
    local v, one, i, j;
	if IsBound(rep.fixedPoints) then
		return rep.fixedPoints;
	fi;
    if IsTrivial(rep.group) then
        return(IdentityMat(rep.dimension,rep.field));
    fi;
    v:=[];
    one:=One(rep.field);
    for i in [1..rep.dimension] do 
        v[i] := Concatenation(List(rep.genimages, m->m[i]));
        for j in [0..Length(rep.genimages)-1] do 
            v[i][i+rep.dimension*j] := v[i][i+rep.dimension*j] - one;
        od;
    od;
    v:=NullspaceMat(v);
    rep.fixedPoints:=v;
    return(v);
end ;

InstallGlobalFunction( FixedQuotient, function(rep)
    local onemat, v;
    if IsBound(rep.fixedQuotient) then
	    return rep.fixedQuotient;
	fi;
    onemat:=SafeIdentityMat(rep.dimension, rep.field);
    v:=Concatenation(List(rep.genimages, g->g - onemat));
    v:=SafeBaseMat(v);
    rep.fixedQuotient:=v;
    return v;
end );

InstallGlobalFunction( SubFixedQuotient, function(rep,u)
    local onemat, v;
    onemat:=SafeIdentityMat(rep.dimension, rep.field);
    v:=Concatenation(List(rep.genimages, g->SafeMatrixMult(u, g - onemat, rep.dimension)));
    return SafeBaseMat(v);
end );

InstallGlobalFunction( BrauerTraceImage, function(rep, p)
    local resp, resh, maxsub, image, h, M;
    resp := RestrictedRep( rep.group, p, rep );
    maxsub:=MaximalSubgroups(p);
    image:=[];
    for h in maxsub do
        M := RelativeTrace(p,h,resp);
        resh := RestrictedRep( resp.group, h, resp );
        Append( image, SafeMatrixMult( FixedPoints@(resh), M, rep.dimension ) );
    od;
    return( SafeBaseMat( image ) );
end );

InstallGlobalFunction( Spin, function(rep, veclist)
    return rep.operations.Spin( rep, veclist );
end );

InstallGlobalFunction( GroupSpin, function(rep, veclist )
    local basis, oldlength, newvecs, g,v;
    basis:=List(SafeBaseMat(veclist));
    oldlength:=Length(basis)-1;
    while Length(basis) > oldlength do
        oldlength:=Length(basis);
        newvecs:=[];
        for g in rep.genimages do
            for v in basis do
                Add(newvecs, v*g);
            od;
        od;
        Append(basis,newvecs);
        basis:=List(SafeBaseMat(basis));
    od;
    return basis;
end );
