global atom_name
atom_name='Rb';
units_and_constants;

% Here are 6 examples for usage of the code.
%
% These are based substantially on examples and results given by
% J. Pritchard, in his PhD Thesis (Durham University 2009).
%
% O. Firstenberg 2012, Harvard University HQOC / MIT

switch 6 % example number (1-6)
    
        
%% Example I: Energy levels
case 1
    figure;
    E_ref=0; unit_factor=(1/hbar/THz);
    N=(5:80);
    L=0; J=1/2;        E_binding=RyM./(nstar(N,L,J).^2);    plot(N,(E_ref-E_binding)*unit_factor,'o-');hold all;
    L=1; J=1/2;        E_binding=RyM./(nstar(N,L,J).^2);    plot(N,(E_ref-E_binding)*unit_factor,'o-');
    L=1; J=3/2;        E_binding=RyM./(nstar(N,L,J).^2);    plot(N,(E_ref-E_binding)*unit_factor,'o-');
    L=2; J=3/2;        E_binding=RyM./(nstar(N,L,J).^2);    plot(N,(E_ref-E_binding)*unit_factor,'o-');
    L=2; J=5/2;        E_binding=RyM./(nstar(N,L,J).^2);    plot(N,(E_ref-E_binding)*unit_factor,'o-');
    L=3; J=5/2;        E_binding=RyM./(nstar(N,L,J).^2);    plot(N,(E_ref-E_binding)*unit_factor,'o-');
    L=3; J=7/2;        E_binding=RyM./(nstar(N,L,J).^2);    plot(N,(E_ref-E_binding)*unit_factor,'o-');
    set(gca,'yscale','log');ylabel('Transission frequency [THz]'); xlabel('n');title('Level structure');
    legend('nS(1/2)','nP(1/2)','nP(3/2)','nD(3/2)','nD(5/2)','nF(5/2)','nF(7/2)',2);

%% Example II: Dipole lengths (radial matrix elements) from 5P3/2 state
case 2
    figure;
    for L2=[0 2]
        N2=30:10:100;     J2=L2+1/2;
        N1=N2*0+5; L1=1;  J1=3/2;
        R=Rdipole_table(N1,L1,J1,N2,L2,J2);
        xplot=nstar(N2,L2,J2);yplot=abs(R);
        plot(log10(xplot),log10(yplot),'o-');I=round(linspace(1,length(xplot),5));set(gca,'XTick',log10(xplot(I)),'XTickLabel',round(xplot(I)),'YTick',log10(sort(yplot(I))),'YTickLabel',round(100*sort(yplot(I)))/100);
        hold all;
    end
    grid on;xlabel('n^*');ylabel('Radial dipole matrix element [a.u.]'); title('Radial matrix elements'); legend('5P_{3/2}->nS_{1/2}','5P_{3/2}->nD_{5/2}');
    
%% Example III: Scalar polarizability
case 3
    figure;colors=jet(3);
    for ind=1:3
        switch ind
            case 1, L1=0; J1=1/2; M1=1/2;
            case 2, L1=2; J1=3/2; M1=1/2;
            case 3, L1=3; J1=5/2; M1=1/2;
        end
        N1_vec=(10:10:80);
        alpha0s=0*N1_vec;
        for ind_N1=1:length(N1_vec);
            N1=N1_vec(ind_N1);
            E0=nstar(N1,L1,J1).^(-2);
            N2s_distance=5;
            N2s=(N1-N2s_distance:N1+N2s_distance); L2s=[L1-1,L1-1,L1+1,L1+1]; J2s=[L1-1.5,L1-0.5,L1+0.5,L1+1.5];
            % The "rulses" are: (i) dJ=0,+1,-1   (ii) J2 not negative, and can support M2=M1, (iii) L2 not negative
            Irulse=(abs(J2s-J1)>1 | J2s<abs(M1) | L2s<0);L2s(Irulse)=[];J2s(Irulse)=[];
            [kuku,L2s]=meshgrid(N2s,L2s);L2s=L2s(:);[N2s,J2s]=meshgrid(N2s,J2s);N2s=N2s(:);J2s=J2s(:);
            IlargeL=(    L2s>=N2s);N2s(IlargeL)=[];L2s(IlargeL)=[];J2s(IlargeL)=[];
            reduced_sqr=abs(reduced_LJM(L1,J1,M1,L2s,J2s,M1)).^2;
            N1s=N2s*0+N1;
            Rdipole12=Rdipole_table(N1s,L1,J1,N2s,L2s,J2s);
            alpha0=sum(reduced_sqr.*abs(Rdipole12).^2./(E0-nstar(N2s,L2s,J2s).^(-2)));
            alpha0=2*alpha0*(a0*qe)^2/(RyM);   % factor of 2 due to definition: alpha0*E^2/2.
            alpha0=alpha0/hbar/MHz*(V^2/cm^2); % convert to MHz*v^-2*cm^2;
            alpha0s(ind_N1)=alpha0;
        end        
        xplot=nstar(N1_vec,L1,J1);yplot=alpha0s;
        plot(log10(xplot),yplot./(xplot.^6),'o');set(gca,'XTickLabel',round(10.^get(gca,'XTick')));hold all;
        plot(log10(xplot),scalar_pol(N1_vec,L1,J1,1/2)./(xplot.^6),'.-');
        xlabel('n^*');ylabel('\alpha_0/(n^*^6) [MHz V^{-2} cm^2]');hold all;
    end
    legend('S(1/2),m=1/2','Analytic fit','D(3/2),m=1/2','Analytic fit','F(5/2),m=1/2','Analytic fit');
    title('Scalar polarizability');
    
%% Example VI: VdV interaction -- C6 vs. n
case 4
    figure;colors=jet(4); 
    geom.type='free space';geom.angle=0;
    for great_ind=1:4
        switch great_ind
            case 1, L1=0*[1 1]; J1=1/2*[1 1];  M1=1/2*[1 1]; 
            case 2, L1=2*[1 1]; J1=3/2*[1 1];  M1=1/2*[1 1];           
            case 3, L1=2*[1 1]; J1=3/2*[1 1];  M1=3/2*[1 1];
            case 4, L1=2*[1 1]; J1=5/2*[1 1];  M1=3/2*[1 1];
        end        
        N1_vec=[23:1:83];
        VdV=N1_vec*0;
        for ind=1:length(N1_vec)
            N1=N1_vec(ind)*[1 1]; 
            VdV(ind)=pair_interaction_old(N1,L1,J1,M1,geom,[]);
        end        
        subplot(1,3,1);        xplot=N1_vec;yplot=VdV./(nstar(N1_vec,L1(1),J1(1)).^11);
        plot(xplot,yplot,'o-','color',colors(great_ind,:),'linewidth',2);hold on;xlabel('n');ylabel('VdV [a.u.], scaled with n*^{11}');
        subplot(1,3,2);        xplot=N1_vec;yplot=-VdV;
        plot(xplot,yplot,'o-','color',colors(great_ind,:),'linewidth',2);hold on;xlabel('n');ylabel('C_6 [a.u.]');
        subplot(1,3,3);        xplot=N1_vec;yplot=-VdV*2*Ry*a0^6/hbar/(GHz*um^6);
        plot(xplot,yplot,'o-','color',colors(great_ind,:),'linewidth',2);hold on;xlabel('n');ylabel('C_6 [GHz \mum^6]');
    end    
    legend('S(1/2),m=1/2','D(3/2),m=1/2','D(3/2),m=3/2','D(5/2),m=3/2',2);

%% Example V: VdV interaction -- C6 vs. theta
case 5
    figure;colors=jet(3);
    N1=46*[1 1]; L1=2*[1 1]; J1=5/2*[1 1];
    geom.type='free space';th_vec=[0:0.1:pi];
    for great_ind=1:3
        switch great_ind
            case 1,  M1=1/2*[1 1];
            case 2,  M1=3/2*[1 1];
            case 3,  M1=5/2*[1 1];
        end        
        VdV=th_vec*0;
        for ind=1:length(th_vec)
            geom.angle=th_vec(ind);
            VdV(ind)=pair_interaction_old(N1,L1,J1,M1,geom,[]);
        end
        xplot=th_vec;yplot=VdV./(nstar(N1(1),L1(1),J1(1)).^11);
        plot(xplot,yplot,'.-','color',colors(great_ind,:),'linewidth',2);hold on;xlabel('\theta');ylabel('VdV shift [a.u.] (scaled with n*^11)');
    end
    legend('(46)D(5/2),m=1/2','(46)D(5/2),m=3/2','(46)D(5/2),m=5/2',2);
    
%% Example VI: Pair interation vs. R
case 6
    figure;colors=jet(4);
    L1=[0 0]; J1=1/2*[1 1];  M1=1/2*[1 1];
    geom.type='free space';geom.angle=0;
    for great_ind=1:4
        switch great_ind
            case 1, N1=30*[1 1];
            case 2, N1=40*[1 1];
            case 3, N1=60*[1 1];
            case 4, N1=80*[1 1];
        end
        R_vec=logspace(log10(0.5),1,40)*um;
        dW=R_vec*0;
        for ind=1:length(R_vec)
            R=R_vec(ind);
            dW(:,ind)=pair_interaction(N1,L1,J1,M1,geom,R/a0);
        end
        xplot=R_vec/um;
        yplot=abs(dW)*2*Ry/hbar/Hz;  % 2*Ry to convert from a.u. to Hz
        plot(xplot,yplot,'.-','color',colors(great_ind,:));
        hold on
        xlabel('R[\mum]');ylabel('|dW| [Hz]');
        set(gca,'yscale','log','xscale','log');
    end
    legend('(30)S(1/2),m=1/2','(40)S(1/2),m=1/2','(60)S(1/2),m=1/2','(80)S(1/2),m=1/2');
end

%% References
% [1] Comparat & Pillet, JOSA B, 27 A208 (2010).
% [2] Li, Mourachko, Noel, Gallagher, PRA 67, 052502 (2003).
% [3] Han, ..., Gallagher, PRA 74, 054502 (2006).
% [4] Devoret, Girvin, Schoelkopf, Ann. Phy. (Leipzig) 16, 767-779 (2007).
% [5] Kaulakys, JPhysB, 28, 4963-4971 (1995)
% [6] Beterov, ..., Entin, PRA 79, 052504 (2009).
% [7] Kamta,... Oumarou, JPhysB 31 (1998) 963–997.
% [8] O'Sullivan & Stoiche?, PRA 31, 2718 (1985).
% [9] Reinhard, ..., Raithel PRA 75, 032712 (2007).