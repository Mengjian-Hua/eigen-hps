function f = mtimes(f, c)
%*   Scale a SURFACEFUN.
%   c*F or F*c multiplies a SURFACEFUN F by a scalar c.
%
%   See also TIMES.

if ( ~isa(f, 'surfacefun') )
    % Ensure F is the SURFACEFUN:
    f = mtimes(c, f);
    return
elseif ( isa(c, 'surfacefun' ) )
    % MTIMES should not be used to multiply two SURFACEFUNs:
    error('SURFACEFUN:mtimes:twofuns', ...
        ['Cannot multiply two surfacefuns with ''*''. ', ...
         'Did you mean ''.*''?\n'])
elseif ( isnumeric(c) && isscalar(c) )
    % Multiply SURFACEFUN F by scalar c:
    n = size(f, 2);
    for i = 1:n
        f(:,i).vals = cellfun(@(vals) c*vals, f(:,i).vals, 'UniformOutput', false);
    end
elseif (isnumeric(c) && isa(f,'surfacefun'))
    if(all(size(c) > 1))
        m = size(c, 1);
        n = size(f, 2);
        f_temp = f*0;
        for i = 1:m
            for j = 1:n
                f_temp(:,i) = f_temp(:,i) + c(j,i)*f(:,i);
            end
        end
        f = f_temp;
    else 
        n = size(f, 2);
        f_temp = f(:,1)*0;
            for j = 1:n
                f_temp = f_temp + c(j)*f(:,j);
            end
        f = f_temp;
    end
else
    error('SURFACEFUN:mtimes:invalid', 'c must be a scalar.')
end

end
