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
R_vec = linspace(2,15,5)*um;

% grid setup
dim = 11;
particles = {[-1 1], 'p'; [1 1], 'p'; [-1 -1], 'p';  [1 -1], 'p';  [0 0], 'a'; [3 1], 'p'; [3 -1], 'p'; [2 0], 'a'; [1 -3], 'p';  [3 -3], 'p'; [2 -2], 'a';};
n_list = zeros(dim, 1);
l_list = zeros(dim, 1);
combinations = (dec2bin(2^dim-1:-1:0)-'0')*2-1;

for r_ind = 1:length(R_vec)
    for np_ind = 30:4:140
        for na_ind = 30:4:140
            Q = zeros(dim);
            for p_ind=1:dim
                if (particles{p_ind,2} == 'p')
                    n_list(p_ind) = np_ind;
                end
                if (particles{p_ind,2} == 'a')
                    n_list(p_ind) = na_ind;
                    l_list(p_ind) = 1;
                end
            end
            for p1=1:dim
                for p2=1:dim
                    if p1~=p2
                        dist = norm(particles{p1, 1} - particles{p2, 1});
                        Q(p1, p2) = eConv*pair_interaction([n_list(p1), n_list(p2)], [l_list(p1), l_list(p2)], [1/2, 1/2], [1/2, 1/2], geom, dist*R_vec(r_ind)./a0);
                    end
                end
            end
            energies = zeros(2^dim, 1);
            for assign=1:2^dim
                binary = combinations(assign,:);
                energies(assign) = binary*(Q*transpose(binary));
            end
            [sorted_energies,sorted_idx]= sort(energies);
            sorted_combinations = combinations(sorted_idx,:);
            
            shouldbe = 0;
            for assign=1:2^dim
                config = sorted_combinations(assign,:);
                parity = sum(config([1 2 3 4]));
                parity1 = sum(config([2 4 6 7]));
                parity2 = sum(config([4 7 9 10]));
                if (abs(parity) == 2 && abs(parity1) == 2 && abs(parity2) == 2)
                    shouldbe = shouldbe + 1;
                else
                    break;
                end
            end
            
            if shouldbe >= 8
                fprintf("%f", shouldbe);
                %config = sorted_combinations(9,:);
                %parity = sum(config([1 2 4 7]));
                %fprintf("%f",abs(parity));
                ratio = pair_interaction([n_list(3), n_list(1)], [l_list(3), l_list(1)], [1/2, 1/2], [1/2, 1/2], geom, R_vec(r_ind)/a0)./pair_interaction([n_list(1), n_list(1)], [l_list(1), l_list(1)], [1/2, 1/2], [1/2, 1/2], geom, R_vec(r_ind)*sqrt(2)/a0);
                a_strength_ratio = pair_interaction([n_list(3), n_list(3)], [l_list(3), l_list(3)], [1/2, 1/2], [1/2, 1/2], geom, R_vec(r_ind)*sqrt(2)/a0)./pair_interaction([n_list(3), n_list(1)], [l_list(3), l_list(1)], [1/2, 1/2], [1/2, 1/2], geom, R_vec(r_ind)/a0);
                p_cross_talk = pair_interaction([n_list(1), n_list(1)], [l_list(1), l_list(1)], [1/2, 1/2], [1/2, 1/2], geom, R_vec(r_ind)*sqrt(2)*2/a0)./pair_interaction([n_list(1), n_list(1)], [l_list(1), l_list(1)], [1/2, 1/2], [1/2, 1/2], geom, R_vec(r_ind)*sqrt(2)/a0);
                fprintf("(X: %f, Delta: %f, a-a interaction: %f, Cross_talk: %f, P: %f, A: %f, R: %0.20f),\n", ratio, (sorted_energies(9) - sorted_energies(8)), abs(a_strength_ratio), abs(p_cross_talk), np_ind, na_ind, R_vec(r_ind));
            end
        end
    end
end