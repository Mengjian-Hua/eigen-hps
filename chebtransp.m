function [an]=chebtransp(f,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  -1-----------------------1    n points, including both endpoints.
%  f1                       fn
%
%  Compute Chebyshev coefficients from function values given 
%  at "practical" nodes xn on [-1,1]. 
%
%  Suppose n = 5: data is [f1, f2, f3, f4, f5]
%  where fj = f(xn(j)):
%
%  We want to unfold this as an even extension on an interval
%  of twice the size (in "theta" coordinates extending from 
%  [0,pi] to [0,2 pi]).
%
%  A little thought shows that the layout should be:
%
%    [f5 f4 f3 f2 f1 f2 f3 f4 ] where we DON'T include 
%
%  the right end-point (2 pi) in the discretization.
%
%  Why is it in the above ordering? Because the map from [0,pi]
%  to [-1, 1] under the map "cosine" identifies cos(0) with the
%  RIGHT endpoint 1 and cos(pi) with the left endpoint -1.
%
%  After laying out the code as above, giving the even extension,
%  it's simply a metter of computing the Fourier transform, scaling
%  by the quadrature weight 1/(n-1) and dividing the zero an(1) by 2
%  so that we are consistent with the convention
%
%    f(x) = c0/2 + sum_{n=1}^\infty  c_n T_n(x)
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n2 = 2*(n-1);
work = 0*ones(n2,1);
work(1) = f(n);
work(n) = f(1);
for k = 2:n-1
work(k)=f(n-k+1);
work(n2-k+2)=work(k);
end
an = 0*f;
ww = fft(work);
scalf = 1.0d0/(n-1);
for k = 1:n
an(k)=ww(k)*scalf;
end
an(1)=an(1)/2;
