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
R_vec = linspace(2, 2, 1)*um;

% grid setup
dim = 8;
particles = {[0 0],'a';[1 0],'p';[0 1],'p';[-1 0],'p';[0 -1],'p';[-1 -1],'a';[-2 -1],'p';[-1 -2],'p';[1 -1],'a';[2 -1],'p';[1 -2],'p'};
%particles = {[0 0],'p';[2 0],'p';[0 2],'p';[2 2],'p';[1 1],'a'};

n_list = zeros(dim, 1);
l_list = zeros(dim, 1);
combinations = transpose((dec2bin(2^(dim)-1:-1:0)-'0')*2-1);

for r_ind = 1:length(R_vec)
    for np_ind = 102:102
        for na_ind = 86:86
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
                        if(particles{p1,2} == 'a' && particles{p2,2} == 'p')
                            Q(p1, p2) = 2;
                        end
                        if(particles{p2,2} == 'a' && particles{p1,2} == 'p')
                            Q(p1, p2) = 2;
                        end
                        if(particles{p2,2} == 'p' && particles{p1,2} == 'p')
                            Q(p1, p2) = 1;
                        end
                            
    
                    end
                end
            end
            Q(2, 7) = 0;
            Q(3, 7) = 0;
            Q(2, 8) = 0;
            Q(3, 8) = 0;
            Q(7, 2) = 0;
            Q(7, 3) = 0;
            Q(8, 2) = 0;
            Q(8, 3) = 0;
            Q(1, 7) = 0;
            Q(1, 8) = 0;
            Q(7,1) = 0;
            Q(8,1) = 0;
            Q(6, 2) = 0;
            Q(6, 3) = 0;
            Q(2,6) = 0;
            Q(3,6) = 0;
            
            energies = zeros(length(combinations), 1);
            for assign=1:length(combinations)
                binary = combinations(:,assign);
                for prime1=1:dim
                    for prime2=1:dim
                        energies(assign) = energies(assign) + Q(prime1, prime2)*binary(prime1)*binary(prime2);
                    end
                end%energies(assign) = transpose(binary)*(Q*binary);
            end
            [sorted_energies,sorted_idx]= sort(energies);
            sorted_combinations = combinations(:,sorted_idx);
            shouldbe8 = 0;
            for counter=1:32
                config = sorted_combinations(:,counter);
                parity = sum(config([2 3 4 5]));
                parity1 = sum(config([4 5 7 8]));
                fprintf('%f %f %f %f \n', counter, parity, parity1, sorted_energies(counter));
            %if (abs(parity) == 2)
                    %shouldbe8 = shouldbe8 + 1;
                %else
                    %break;
                %end
            end
        end
    end
end