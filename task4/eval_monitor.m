function f = eval_monitor(lambda,pdo,dom)
    
    if(isfield(pdo,'c'))
        pdo.b = pdo.b-lambda;
    else
        pdo.b = -lambda;
    end
    
    L = surfaceop(dom, pdo, 0);
    % L.rankdef = true;
    build(L)
    % D2N1 = L.patches{1,1}.child1.D2N;
    % D2N2 = L.patches{1,1}.child2.D2N;
    f = eigs(L.patches{1}.A,1,"smallestabs");

end