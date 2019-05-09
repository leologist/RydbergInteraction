function str = StateLabelString(state)
%str = StateLabelString(state)
%   state = [n,L,J,m];
%   

n = state(1);
L = state(2);
J = state(3);
m = state(4);

Lstrs = {'S','P','D','F'};

str = sprintf('%d%s_{%d/2,%d/2}', n, Lstrs{L+1}, J*2, m*2);
end

