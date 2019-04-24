function y=Rdipole_table(n1,L1,J1,n2,L2,J2)
% Yields the radial dipole matrix element between states
% |n1 L1 J1> and |n2 L2 J2> in units of a.u.
%
% Vectorized (the input arrays should share the same size, or be scalars).
%
% Uses a pre-calculated table of values (.mat file)
%
% IMPORTANT: The first time the function runs in each session, it takes a
% few minutes to "hash" the data in the table, during which matlab will
% become unresponsive.
%
% In the absence of the table file, it can be constructed 
% by running the funciton with n1=999. For "safety" reasons, the "save"
% command is currently commented out. Remeber to uncomment it, or save
% manually
%
% O. Firstenberg, Harvard University HQOC ; MIT, 2012.
%

global Rdipole_global_table
global qn_Rdipole_global_table

nsMax = 200; LsMax = 4; 
if isempty(qn_Rdipole_global_table)
    if n1~=-999
        disp('Loading global table for Rdipole (once in a MATLAB session).');
        load(sprintf('Rdipole_table_n%d_n%d_L%d_L%d_J2_J2',nsMax,nsMax,LsMax,LsMax),'Rdipole_global_table', 'qn_Rdipole_global_table');
    else
        % Rebuild the table
        ns=1:nsMax; Ls=0:LsMax;  Js=[1,-1]/2;
        [ns,Ls,Js]=meshgrid(ns,Ls,Js);
        ns=ns(:);Ls=Ls(:);Js=Js(:)+Ls;
        Inegative=(Ls>=ns | Js<0); ns(Inegative)=[];Ls(Inegative)=[];Js(Inegative)=[];
        [n1s,n2s]=meshgrid(ns);[L1s,L2s]=meshgrid(Ls);[J1s,J2s]=meshgrid(Js);
        Itriu=triu(n1s*0+1,1);Itriu=(Itriu(:)==1);
        qn=[n1s(Itriu) L1s(Itriu) J1s(Itriu) n2s(Itriu) L2s(Itriu) J2s(Itriu)];
        Irules=( (abs(J2s(Itriu)-J1s(Itriu))>1) | (abs(L2s(Itriu)-L1s(Itriu))>1) );
        qn(Irules,:)=[];
        clear Irules Itriu Inegative ns n1s n2s Ls L1s L2s Js J1s J2s
    
        % Find each radial dipole matrix element
        Rdipole_global_table=qn(:,1)*0;
        parts=[1 round(linspace(2,length(Rdipole_global_table),99)) length(Rdipole_global_table)+1];        
        for ind=1:100
            disp([num2str(ind) '%']);
            I=parts(ind):(parts(ind+1)-1);
            if ~isempty(I), Rdipole_global_table(I)=Rdipole(qn(I,1),qn(I,2),qn(I,3),qn(I,4),qn(I,5),qn(I,6)); end
        end
        
        disp('Hashing table. MATLAB will become unresponsive for a while (wait for message)...');
        hash_ind=qn(:,1)*1e7+qn(:,4)*1e4+qn(:,2)*1e3+qn(:,5)*1e2+(qn(:,3)+1/2)*1e1+(qn(:,6)+1/2)*1e0;
        qn_Rdipole_global_table=sparse(max(hash_ind),1);
        qn_Rdipole_global_table(hash_ind)=(1:length(hash_ind))';
        clear qn hash_ind;
        disp('Done hashing.')
        % remember to save: 
        save(sprintf('Rdipole_table_n%d_n%d_L%d_L%d_J2_J2',nsMax,nsMax,LsMax,LsMax),'Rdipole_global_table', 'qn_Rdipole_global_table');
    end

end
hash_ind1=n1*1e7+n2*1e4+L1*1e3+L2*1e2+(J1+1/2)*1e1+(J2+1/2)*1e0;
hash_ind2=n2*1e7+n1*1e4+L2*1e3+L1*1e2+(J2+1/2)*1e1+(J1+1/2)*1e0;
table_ind=qn_Rdipole_global_table(hash_ind1);
Izeros=(table_ind==0);
table_ind(Izeros)=qn_Rdipole_global_table(hash_ind2(Izeros));
y=Rdipole_global_table(table_ind);
y(Izeros)=conj(y(Izeros));
if ~all(size(y)==size(hash_ind1)), y=y.'; end
