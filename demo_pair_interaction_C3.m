global atom_name
atom_name = 'Rb';
units_and_constants;
geom.type='free space'; geom.angle=0;

figind = 12; % figure index

% n1 = 100; l1 = 0; j1 = 1/2; m1 = 1/2;
n1 = 69; l1 = 1; j1 = 3/2; m1 = 3/2;

n2 = 70; l2 = 0; j2 = 1/2; m2 = 1/2;
% n2 = n1; l2 = l1; j2 = j1; m2 = m1;

flipflop = any([n1-n2, l1-l2, j1-j2, m1-m2]);


R_vec = logspace(log10(2.5),4,100)*um;

eConv = 2*Ry/hbar/Hz; % convert from energy in atomic units to SI freq

Hamsizes = zeros(size(R_vec));
energies_relevant = [];
overlaps_relevant = [];
minEs = [];

tic;
for ind = 1:length(R_vec)
    [~, eigval, ov, info] = ...
        pair_interaction([n1, n2], [l1, l2], [j1, j2], [m1, m2], geom, R_vec(ind)/a0, 9,1e-6);
    [temp, idx] = max(abs(ov).^2);
    energies_relevant = [energies_relevant; reshape(eigval(idx)*eConv, size(idx))];
    temp2 = eigval(idx);
    meanE = mean(temp2);
    [~, Imin] = min(abs(temp2-meanE));
    minEs = [minEs; temp2(Imin)*eConv];
    overlaps_relevant = [overlaps_relevant; temp];
    Hamsizes(ind) = length(eigval);
end
toc;

%% printing out asymptotic wavefunctions

numAsymptotes = size(overlaps_relevant, 2);
legend_labels = cell(1, numAsymptotes);

fprintf('======\n')
fprintf('Asymptotic Wavefunctions\n')
fprintf('% 32s: ','State');
fprintf('    % 5i', 1:numAsymptotes);
fprintf('\n');
for ind = 1:numAsymptotes
    legend_labels{ind} = ['State ' num2str(ind)];
    fprintf('% 15s, % 15s: ', StateLabelString(info.state_labels(ind,1:2:end)),...
        StateLabelString(info.state_labels(ind,2:2:end)));
    fprintf('    % +4.2f', full(info.overlap_vectors(ind, :)));
    fprintf('\n');
end


%% Plotting results


C3_calc = min(abs(eig(info.HamResSect)))*a0^3*eConv;

figure(figind);

subplot(2,2,1)
plot(R_vec/um, abs(energies_relevant),'.-');
hold on
plot(R_vec/um, abs(C3_calc./(R_vec).^3),'--m');
% plot(R_vec/um, abs(C6nsns(n1))./R_vec.^6,'--r');
% plot(R_vec/um, abs(C5npnp(n1)./R_vec.^5),'--r');
% plot(R_vec/um, abs(C6npnp(n1)./R_vec.^6),'--c');
% plot(R_vec/um, abs(C5npnp(n1)./R_vec.^5 + C6npnp(n1)./R_vec.^6),'--r');
plot(R_vec/um, abs(minEs),'--g');
hold off
set(gca,'xscale','log','yscale','log');
xlabel('R (\mum)');
ylabel('|Energy of Interaction| (Hz)');
legend(legend_labels)
grid on



y = mean(energies_relevant,2)*ones(1,numAsymptotes);

subplot(2,2,2);
plot(R_vec/um, abs(energies_relevant-y),'.-',R_vec/um, abs(y(:, 1)));
hold on
plot(R_vec/um, abs(C3_calc./(R_vec).^3),'--m');
hold off
set(gca,'xscale','log','yscale','log');
legend([legend_labels, 'average', 'C_3/R^3 asymp']);
xlabel('R (\mum)');
ylabel('|Energy of Interaction| (Hz)');
title('|Differences from average|');
grid on

subplot(2,2,3);
plot(R_vec/um, overlaps_relevant,'.-');
set(gca,'xscale','log','yscale','lin','ylim',[0,1.1]);
xlabel('R (\mum)');
ylabel('|overlap|^2');
title('eigenstate(R) overlap with eigenstate(R=\infty)');
legend(legend_labels,'location','southeast');

subplot(2,2,4)
plot(R_vec/um, Hamsizes,'.-');
set(gca,'xscale','log')
xlabel('R (\mum)');
ylabel('dimension of Hilbert space diagonalized');
