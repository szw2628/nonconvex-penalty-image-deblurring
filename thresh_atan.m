function x = thresh_atan(y, lambda, a)
% x = thresh_atan(y, lambda, a)
%
% Threshold function induced by the arctangent penalty
%
%   x = arg min_x { 0.5 * (y - x)^2 + lambda * phi(x, a) }
%   where
%   phi(x, a) =  2/(a*sqrt(3)) * (atan((2*a*abs(x)+1)/sqrt(3)) - pi/6)
%   and
%   0 <= a <= 1/lambda
%
% INPUT
%   y : data (scalar or multidimensional array)
%   lambda : threshold (scalar or multidimensional array)
%   a : convexity parameter (a >= 0)
%       If a = 0 then it reduces to the soft-threshold function.
% OUTPUT
%   x : output of thresholding


if (a > 1/lambda)
    disp('a should be less than 1/lamba.')
    x = [];
    return
end
if (a < 0)
    disp('a should be non-negative.')
    x = [];
    return
end

if ( a < 1e-10 )
    
    x = soft(y, lambda);
    
else
    b = 1 - a.*abs(y);
    
    i = find( b == 0 );
    
    c1 = b.^3./(27.*a^3) - (b.^2)./(6.*a^3) - (abs(y) - lambda)./(2.*a^2) ;
    
    c2 = b.^2./(9.*a^2) - b./(3.*a^2);
    
    c3 = ((c1.^2 - c2.^3).^(1/2) - c1).^(1/3);
    
    z = c3 - b./(3.*a) + c2./c3;
    
    if isempty(i) ~= 1
        z(i) = ( (abs(y(i)) - lambda )./a^2 ) .^ (1/3);
    end
    
    x = abs(z) .* sign(y) .*( (abs(y)-lambda) >= 0 );
    
end

end

function x = soft(y, lambda)
% x = soft(y, lambda)
%
% SOFT THRESHOLDING
% for real or complex data.
%
% INPUT
%   y : data (scalar or multidimensional array)
%   lambda : threshold (scalar or multidimensional array)
%
% OUTPUT
%   x : output of soft thresholding
%
% If x and lambda are both multidimensional, then they must be of the same size.

x = max(1 - lambda./abs(y), 0) .* y;
end
