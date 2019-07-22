global atom_name
atom_name = 'Rb';

MainWorkspaceDir = pwd;
addpath([MainWorkspaceDir, filesep, 'atomic_data']);

units_and_constants;
geom.type='free space'; geom.angle=0;


figind = 11;

% n1 = 60; l1 = 0; j1 = 1/2; m1 = 1/2; %100S_{1/2}
% n1 = 40; l1 = 1; j1 = 3/2; m1 = 3/2; % 100P_{3/2}
n1 = 30; l1 = 2; j1 = 5/2; m1 = 5/2; % 100P_{3/2}
n2 = n1; l2 = l1; j2 = j1; m2 = m1;
eConv = 2*Ry/hbar/Hz; % convert from energy in atomic units to SI freq

R_vec = logspace(log10(2.5),3,60)*um;
Vryd = nan(length(R_vec),1);
Vryd2 = nan(length(R_vec),1);

for ind = 1:length(R_vec)
    tic;
    Vryd(ind,:) = pair_interaction([n1, n2], [l1, l2], [j1, j2], [m1, m2], geom, R_vec(ind)/a0);
    fprintf('%0.2e um, time = %0.2f sec', R_vec(ind)/um, toc);
    fprintf('\n');
%     tic;
%     Vryd2(ind,:) = pair_interaction_old([n1, n2], [l1, l2], [j1, j2], [m1, m2], geom, R_vec(ind)/a0);
%     fprintf(', %0.2f sec\n', toc);
end

%%

% Rb data from Singer et al (doi:10.1088/0953-4075/38/2/021)

C6nsns = @(n) -n.^11.*(1.197e1 - 8.486e-1*n + 3.385e-3*n.^2)*a0^6*eConv;
C5npnp = @(n) -n.^8.*(0.922 + 7.903-2*n  -0.041e-2*n.^2)*a0^5*eConv;
C6npnpfun = @(n, c0, c1, c2) -n.^11.*(c0*1e-1 - c1*1e-1*n + c2*1e-4*n.^2)*a0^6*eConv;
C6npnp = @(n) C6npnpfun(n, 3.575, -0.183, 0.816);
% C6npnp = @(n) C6npnpfun(n, 5.461, -1.133, 5.476);

C6ref = 0;
if n1 == n2 && l1 == l2
    if l1 == 0 
        C6ref = C6nsns(n1);
    elseif l1 == 1
        C6ref = C6npnp(n1);
    end
end
    

% C6 estimated from second-order perturbation theory 
C6estimate = pair_interaction_old([n1, n2], [l1, l2], [j1, j2], [m1, m2], geom,[])*a0^6*eConv;

figure(figind);
plot(R_vec/um, abs(Vryd)*eConv,'.-');
hold on
% plot(R_vec/um, abs(Vryd2)*eConv,'o-');
plot(R_vec/um, abs(C6ref./R_vec.^6),'--r');
plot(R_vec/um, abs(C6estimate./R_vec.^6),':g');
hold off, grid on
set(gca,'yscale','log','xscale','log');
xlabel('R [\mum]');
ylabel('Interaction');
legend('Calculated from ED', 'from literature', 'estimated');

