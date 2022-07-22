function P = cheb_projection(n)
    
    vals_to_coeffs = zeros(n-1,n-1);
    coeffs_to_vals = zeros(n-1,n-1);
    
    for i = 1:n-1
        theta = (i-1)*pi/(n-2);
        for j = 1:n-1
            vals_to_coeffs(i,j) = 2/(n-1)*cos(pi*i*(j-1/2)/(n-1));
            coeffs_to_vals(i,j) = cos(j*theta);
        end
    end
    
    P = coeffs_to_vals*vals_to_coeffs;
    
end