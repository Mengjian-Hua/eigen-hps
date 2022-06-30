function f = uminus(f)
%-   Unary minus for a SURFACEFUN.
%   -F negates the SURFACEFUN F.
%
%   G = uminus(F) is called for the syntax '-F'.
%
%   See also UPLUS.
[~,m] = size(f);
    if(m == 1)
        f.vals = cellfun(@uminus, f.vals, 'UniformOutput', false);
    else 
        for i = 1:m
            f(:,i).vals = cellfun(@uminus, f(:,i).vals, 'UniformOutput', false);
        end
    end

end
