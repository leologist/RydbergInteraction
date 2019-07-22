% atom = Rubidium 
global atom_name
atom_name = 'Rb';

% import units and constants 
MainWorkspaceDir = pwd;
addpath([MainWorkspaceDir, filesep, 'atomic_data']);
addpath([MainWorkspaceDir, filesep, 'math']);
addpath([MainWorkspaceDir, filesep, 'util']);
units_and_constants;
geom.type='free space'; geom.angle=0;
eConv = 2*Ry/hbar/Hz;

% distances in m
R_vec = sqrt(8)*um;

a = struct
x = struct
y = struct
c = struct

a.n = 86; a.l = 1; a.j = 1/2; a.m = 1/2;
x.n = 102; x.l = 0; x.j = 1/2; x.m = 1/2;
y.n = 119; y.l = 1; y.j = 3/2; y.m = 1/2;
c.n = 100; c.l = 2; c.j = 3/2; c.m = 1/2;

particles = {[0 0],a;[1 0],x;[0 1],x;[-1 0],x;[0 -1],y;[-1 -1],a;[-2 -1],x;[-1 -2],c;[1 -1],a;[2 -1],x;[1 -2],c};
dim = 11;
Q = zeros(dim);

for p1=1:dim
    for p2=1:dim
        if p1~=p2
            dist = norm(particles{p1, 1} - particles{p2, 1});
            alpha = particles{p1,2};
            beta = particles{p2,2};
            Q(p1, p2) = eConv*pair_interaction([alpha.n, beta.n], [alpha.l, beta.l], [alpha.j, beta.j], [alpha.m, beta.m], geom, dist*R_vec/a0);
        end
    end
end