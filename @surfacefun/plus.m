function h = plus(f, g)
%+   Plus for SURFACEFUN.
%   F + G adds the SURFACEFUN F and G. F and G must have the same domains
%   and discretization sizes. F and G may also be scalars.
%
%   See also MINUS.

if ( isnumeric( f ) )
    h = plus(g, f);
    return
elseif ( isnumeric( g ) )
    h = f;
    for k = 1:size(h.vals,1)
        h.vals{k} = h.vals{k} + g;
    end
elseif ( isa(f, 'surfacefun') && isa(g, 'surfacefun') )
    h = f;
    % TODO: Assume on the same grid for now.
    [~,m] = size(f);
    [~,n] = size(g);
    if(m==n)
        for i = 1:m
            h(:,i).vals = cellfun(@plus, f(:,i).vals , g(:,i).vals, 'UniformOutput', false);
        end
    elseif(m>1 && n==1)
        for i = 1:n
            h(:,i).vals = cellfun(@plus, f(:,i).vals , g.vals, 'UniformOutput', false);
        end
    elseif(m == 1 && n>1)
        h = plus(g,f);
    else
        error('SURFACEFUN:mtimes:invalid', 'dimensions of f and g are not compatible.')
    end
%     h.vals = cellfun(@plus, f.vals , g.vals, 'UniformOutput', false);
end

end
