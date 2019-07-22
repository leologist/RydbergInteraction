global atom_name
atom_name = 'Rb';

MainWorkspaceDir = pwd;
addpath([MainWorkspaceDir, filesep, 'atomic_data']);
addpath([MainWorkspaceDir, filesep, 'math']);
addpath([MainWorkspaceDir, filesep, 'util']);

units_and_constants;
geom.type='free space'; geom.angle=0;

R_vec = 4.24*um;

particles = {[0 2], 'p'; [-1 1], 'p'; [0 1], 'a'; [1 1], 'p'; [-2 0], 'p'; [-1 0], 'a'; [0 0], 'p'; [1 0], 'a'; [2 0], 'p'; [-1 -1], 'p'; [1 -1], 'p'};
Q = zeros(11);
np = 80;
lp = 0;
jp = 1/2;
mp = 1/2;

na = 76;
la = 0;
ja = 1/2;
ma = 1/2;


nlist = zeros(11, 1);
llist = zeros(11, 1);
jlist = zeros(11, 1);
mlist = zeros(11, 1);

for p=1:11
    if (particles{p,2} == 'p')
        nlist(p) = np;
        llist(p) = lp;
        jlist(p) = jp;
        mlist(p) = mp;
    end
    if (particles{p,2} == 'a')
        nlist(p) = na;
        llist(p) = la;
        jlist(p) = ja;
        mlist(p) = ma;
    end
end

H = zeros(11, 1);
H(1) = -5;
%H(2) = 3;
%H(4) = -5.5;
%H(5) = -5.3;
%H(7) = 81;
%H(9) = -5;
%H(10) = 500;
%H(11) = 500;

for p1=1:11
    for p2=p1+1:11
        dist = norm(particles{p1, 1} - particles{p2, 1});
        %fprintf("%0.2f \n", dist);
        Q(p1, p2) = pair_interaction([nlist(p1), nlist(p2)], [llist(p1), llist(p2)], [jlist(p1), jlist(p2)], [mlist(p1), mlist(p2)], geom, R_vec*dist/a0);
        Q(p2, p1) = Q(p1, p2);
    end
end

combinations = (dec2bin(2^11-1:-1:0)-'0')*2-1;

energies = diag(combinations*(Q*transpose(combinations)));

extraterms = transpose(transpose(H)*transpose(combinations));

energies = energies + extraterms;

[energies,idx]= sort(energies);

combinations = combinations([idx],:);

%combinations(:,[1 2 4 7]),2)


