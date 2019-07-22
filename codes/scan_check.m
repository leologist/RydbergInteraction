global atom_name
atom_name = 'Rb';

MainWorkspaceDir = pwd;
addpath([MainWorkspaceDir, filesep, 'atomic_data']);
addpath([MainWorkspaceDir, filesep, 'math']);
addpath([MainWorkspaceDir, filesep, 'util']);

units_and_constants;
geom.type='free space'; geom.angle=0;

R_vec = [10];

eConv = 2*Ry/hbar/Hz;

%dim = 11;
dim = 5;

%particles = {[0 2], 'p'; [-1 1], 'p'; [0 1], 'a'; [1 1], 'p'; [-2 0], 'p'; [-1 0], 'a'; [0 0], 'p'; [1 0], 'a'; [2 0], 'p'; [-1 -1], 'p'; [1 -1], 'p'};
particles = {[0 2], 'p'; [-1 1], 'p'; [0 1], 'a'; [1 1], 'p'; [0 0], 'p'};

spin = 0.5;
l = 0;
momentum = 0.5;

n_list = zeros(dim, 1);
l_list = zeros(dim, 1);
j_list = zeros(dim, 1);
m_list = zeros(dim, 1);

%combinations = [dec2bin(2^9-1:-1:0)-'0' zeros(512,2)];
combinations = (dec2bin(2^5-1:-1:0)-'0')*2-1;

for r_ind = 1:length(R_vec)
    for np_ind = 115:115
        fprintf('Checking physical level %i \n', np_ind);
        for na_ind = 95:95
            %fprintf('%f', na_ind);
            Q = zeros(dim);
            for p_ind=1:dim
                l_list(p_ind) = l;
                m_list(p_ind) = spin;
                j_list(p_ind) = momentum;
                if (particles{p_ind,2} == 'p')
                    n_list(p_ind) = np_ind;
                end
                if (particles{p_ind,2} == 'a')
                    n_list(p_ind) = na_ind;
                    l_list(p_ind) = 1;
                    %j_list(p_ind) = 3/2;
                end
            end
            for p1=1:dim
                for p2=p1+1:dim
                    dist = norm(particles{p1, 1} - particles{p2, 1});
                    Q(p1, p2) = eConv*pair_interaction([n_list(p1), n_list(p2)], [l_list(p1), l_list(p2)], [j_list(p1), j_list(p2)], [m_list(p1), m_list(p2)], geom, R_vec.*dist.*um./a0);
                    Q(p2, p1) = Q(p1, p2);
                end
            end
            
            energies = diag(combinations*(Q*transpose(combinations)));
            [energies,idx]= sort(energies);
            check = 0;
            for ground_states=1:8
                parity = combinations(idx(ground_states),:);
                parity1 = sum(parity([1 2 4 5]));
                %parity2 = sum(parity([2 5 7 10]));
                %parity3 = sum(parity([4 7 9 11]));
                %fprintf("%f %f %f \n", parity1, parity2, parity3);
                if (abs(parity1) == 2)% && abs(parity2) == 2 && abs(parity3) == 2)
                    check = check + 1;
                else
                    break;
                end
            end
            fprintf("%f\n", check);
            %if check == 8  
                %parity = combinations(idx(9),:);
               % parity1 = sum(parity([1 2 4 7]));
                %parity2 = sum(parity([2 5 7 10]));
                %parity3 = sum(parity([4 7 9 11]));
                %fprintf("%f %f %f \n", parity1, parity2, parity3);
                %fprintf("%f \n", parity1);
                %fprintf("%f\n", );
               % strength_ratio = pair_interaction([n_list(3), n_list(3)], [l_list(3), l_list(3)], [j_list(3), j_list(3)], [m_list(3), m_list(3)], geom, 5.*1.414.*um./a0)./pair_interaction([n_list(1), n_list(1)], [l_list(1), l_list(1)], [j_list(1), j_list(1)], [m_list(1), m_list(1)], geom, 5.*1.414*um./a0);
                %fprintf("gap: %.20f, strength_ratio: %f, physical: %f, ancilla: %f\n", energies(9) - energies(8), strength_ratio, np_ind, na_ind);
            %end
            fprintf("\n\n");
        end
    end
end