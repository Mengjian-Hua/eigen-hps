function lambda = root_finder(lambda_0,dom,pdo,tol)
%% This function tries to find the closet eigenvalue near the initial guess
% lambda_0 for the operator pdo defined in the closed domain dom 
% The tolerance for the error in terms of the monitor function is tol.

    % evaluate the monitor function at the initial guess
    f1 = eval_monitor(lambda_0,pdo,dom); 
    
    % find two points close to the initial guess
    diff = 1e-2;
    lambda = lambda_0;
    lambda2 = lambda - diff;
    lambda3 = lambda + diff;
    
    while (f1>tol)
        
        % compute left and right derivatives
        f2 = eval_monitor(lambda2,pdo,dom);
        f3 = eval_monitor(lambda3,pdo,dom);
        
        df1 = (f1 - f2)/(diff);
        df2 = (f3 - f1)/(diff);
            
        % if the left and right derivatives have different signs
        % then lambda, lambda_2, lambda_3 are not on the same piece
        % and we should refine our search. If not, then we do Newton's
        % method until the error is less than the tolerance. 
        
            if(df1*df2<0)
                diff = diff/2;
            else
                lambda = lambda - f1/df1;
                lambda2 = lambda + diff;
                f1 = eval_monitor(lambda,pdo,dom);
            end 
            
            if(f1/df1 < 1e-4)
                break
            end
    end

end

