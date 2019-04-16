%% Workspace script: Establish units and constants variable
%
% O. Firstenberg, Harvard University HQOC ; MIT, 2012.
%%
MainWorkspaceDir = pwd;
addpath([MainWorkspaceDir, filesep, 'atomic_data']);
addpath([MainWorkspaceDir, filesep, 'math']);
addpath([MainWorkspaceDir, filesep, 'util']);

set(0,'DefaultAxesFontSize',16);

%%
global ReducedMassFactor
global atom_name

%% Units
s=1; % second
Hz=2*pi/s;KHz=1000*Hz;MHz=1000*KHz;GHz=1000*MHz;THz=1000*GHz;
kg=1; % kilogram
C=1; % Coulomb
m=1; % meter
cm=1/100;mm=1/1000;um=mm/1000; % centimeter, millimeter, and micrometer

J=kg*m^2/s^2; % Joul
V=J/C; % Volts
F=C/V; % Farad

%% Constants
hbar=1.054571726e-34 *J*s;
epsilon0=8.854187817620e-12 *F/m;
me=9.10938291e-31 *kg;
qe=1.602176565e-19 *C;
c=299792458 *m/s;

k=1/(4*pi*epsilon0);
alpha=k*qe^2/hbar/c; % =1/137
a0=hbar/me/c/alpha; % [m] a0=0.53 Angstrum
Ry=alpha^2/2*me*c^2; % Rydberg constant in [J] == 13.6 [eV]

eV=qe*V; % [Jouls]

%% Atom-specific constants
switch atom_name
    case 'Rb', ma=1.44316060e-25 *kg; % atomic mass for Rb87
    case 'Cs', ma=2.20694650e-25 *kg; % atomic mass for Cs
end
ReducedMassFactor=(1+me/ma);
RyM=Ry/ReducedMassFactor; % Corrected Rydberg energy for the reduced mass