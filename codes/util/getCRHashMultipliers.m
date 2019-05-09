function multi = getCRHashMultipliers(hashprimes)
% multi = getCRHashMultipliers(hashprimes)
%
% get hash multipliers from an array of hash primes
% to use for one-to-one hashing of a list of integers in an algorithm 
% inspired by the Chinese remainder algorithm
% 
% Example:
%   To hash a list of integers (a1, a2, ...)
%   HASH = a1*multi(1) + a2*multi(2) + ...
%   
%   Hash = ai (mod hashprimes(i))
%

N = prod(hashprimes);
Ni = N./hashprimes;

multi = zeros(size(hashprimes));
for ind = 1:length(hashprimes)
    [g, M, ~] = gcd(Ni(ind), hashprimes(ind));
    multi(ind) = M*Ni(ind);
    if g ~= 1
        error('not coprime');
    end
end

end