% Fig8-9.m

% This computer code is from "Interpreting Frequency Responses to Dose-Conserved Pulsatile
% Input Signals in Simple Cell Signaling Motifs", by P. Fletcher, F. Clement, A. Vidal,
% J. Tabak, and R. Bertram, PLoS One, 9:e95613, 2014.

function fig8_9_ab_mismatched()
% This code reproduces the panels for figures 8 and 9 of
%
% Fletcher PA, Clement F, Vidal A, Tabak J, Bertram R (2014)
% Interpreting Frequency Responses to Dose-Conserved Pulsatile Input
% Signals in Simple Cell Signaling Motifs. PLoS ONE 9(4): e95613.
% doi:10.1371/journal.pone.0095613
%
% Figure 8. Input signal dose can affect whether a bell-shaped
% response is observed (panels: figures 1-4)
%
% Figure 9. The kinetics of the linear component affect the peak
% frequency of a bell shaped response (panels: figures 5-8)

% NOTE: There is a minor typo in the paper regarding the b-equation. 
% The b-equation corresponding to the figures would
% be: taub*db/dt=10*F(a)-b.  
% The shapes of the responses are unaffected; this simply changes 
% the vertical scaling of <b>.  See line 221 of this code.

taua=10;

%Signals
T0=100;
f=10.^(-4:0.01:-1); T=1./f;

%hill function
n=2; K=1; Kn=K^n;
F=@(x) x.^n./(x.^n+Kn);

%numerical integral tolerance
tol=1e-6;

%% A0 figure

%Amplitude compensation panels
d0=10; d=d0*ones(size(f));

A0=1; A=A0*T./T0;
[meanAon1, meanb1]=computeMeans();

A0=3; A=A0*T./T0;
[meanAon2, meanb2]=computeMeans();

A0=10; A=A0*T./T0;
[meanAon3, meanb3]=computeMeans();


%plot amplitude compensation panels
figure(1);clf
loglog(f,meanAon1,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
loglog(f,meanAon2,'-','color',[0,0.5,0],'linewidth',2)
loglog(f,meanAon3,'--','color',[0.5,0.5,0],'linewidth',2)
loglog([min(f),max(f)],[1,1],'k:','linewidth',1)
set(gca,'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
set(gca,'ytick',[0.001,0.01,0.1,1,10,100],'yticklabel',{'0.001','0.01','0.1','1','10','100'})
axis([0.0001,0.1,0.1,10])
xlabel('f')
ylabel('<a>_{on,\infty}')
axis square

%plot F(<a>), amp+A0
dF=0.02;


figure(2);clf
semilogx(f,meanb1,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
semilogx(f,meanb2,'-','color',[0,0.5,0],'linewidth',2)
semilogx(f,meanb3,'--','color',[0.5,0.5,0],'linewidth',2)
set(gca,'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
axis([0.0001,0.1,0,5])
xlabel('f');
ylabel('<b>_{\infty}')
axis square

%Duration compensation panels
d0=10; d=d0*T./T0;

A0=1; A=A0*ones(size(f));
[meanAon1, meanb1]=computeMeans();

A0=3; A=A0*ones(size(f));
[meanAon2, meanb2]=computeMeans();

A0=10; A=A0*ones(size(f));
[meanAon3, meanb3]=computeMeans();




%plot amplitude compensation panels
figure(3);clf
loglog(f,meanAon1,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
loglog(f,meanAon2,'-','color',[0,0.5,0],'linewidth',2)
loglog(f,meanAon3,'--','color',[0.5,0.5,0],'linewidth',2)
loglog([min(f),max(f)],[1,1],'k:','linewidth',1)
set(gca,'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
set(gca,'ytick',[0.001,0.01,0.1,1,10,100],'yticklabel',{'0.001','0.01','0.1','1','10','100'})
axis([0.0001,0.1,0.1,10])
xlabel('f')
ylabel('<a>_{on,\infty}')
axis square



figure(4);clf
semilogx(f,meanb1,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
semilogx(f,meanb2,'-','color',[0,0.5,0],'linewidth',2)
semilogx(f,meanb3,'--','color',[0.5,0.5,0],'linewidth',2)
set(gca,'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
axis([0.0001,0.1,0,5])
xlabel('f');
ylabel('<b>_{\infty}')
axis square



%% tau figure

%Amplitude compensation panels
d0=10; d=d0*ones(size(f));
A0=3; A=A0*T./T0;

taua=1;
[meanAon1, meanb1]=computeMeans();

taua=10;
[meanAon2, meanb2]=computeMeans();

taua=100;
[meanAon3, meanb3]=computeMeans();

%plot amplitude compensation panels
figure(5);clf
loglog(f,meanAon1,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
loglog(f,meanAon2,'-','color',[0,0.5,0],'linewidth',2)
loglog(f,meanAon3,'--','color',[0.5,0.5,0],'linewidth',2)
loglog([min(f),max(f)],[1,1],'k:','linewidth',1)
set(gca,'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
set(gca,'ytick',[0.001,0.01,0.1,1,10,100],'yticklabel',{'0.001','0.01','0.1','1','10','100'})
axis([0.0001,0.1,0.1,10])
xlabel('f')
ylabel('<a>_{on,\infty}')
axis square

figure(6);clf
semilogx(f,meanb1,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
semilogx(f,meanb2,'-','color',[0,0.5,0],'linewidth',2)
semilogx(f,meanb3,'--','color',[0.5,0.5,0],'linewidth',2)
set(gca,'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
axis([0.0001,0.1,0,1.5])
xlabel('f');
ylabel('<b>_{\infty}')
axis square

%Duration compensation panels
d0=10; d=d0*T./T0;
A0=3; A=A0*ones(size(f));

taua=1;
[meanAon1, meanb1]=computeMeans();

taua=10;
[meanAon2, meanb2]=computeMeans();

taua=100;
[meanAon3, meanb3]=computeMeans();

%plot amplitude compensation panels
figure(7);clf
loglog(f,meanAon1,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
loglog(f,meanAon2,'-','color',[0,0.5,0],'linewidth',2)
loglog(f,meanAon3,'--','color',[0.5,0.5,0],'linewidth',2)
loglog([min(f),max(f)],[1,1],'k:','linewidth',1)
set(gca,'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
set(gca,'ytick',[0.001,0.01,0.1,1,10,100],'yticklabel',{'0.001','0.01','0.1','1','10','100'})
axis([0.0001,0.1,0.1,10])
xlabel('f')
ylabel('<a>_{on,\infty}')
axis square

figure(8);clf
semilogx(f,meanb1,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
semilogx(f,meanb2,'-','color',[0,0.5,0],'linewidth',2)
semilogx(f,meanb3,'--','color',[0.5,0.5,0],'linewidth',2)
set(gca,'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
axis([0.0001,0.1,0,1.5])
xlabel('f');
ylabel('<b>_{\infty}')
axis square



%%% NESTED HELPER FUNCTIONS
    function [meana1A, meanb]=computeMeans()
        
        alph=A.*(1-exp(-d/taua));
        maxAinf=alph./(1-exp(-T/taua)); %steady state maximum value of "a" at the end of each on-phase
        minAinf=exp(-(T-d)/taua).*maxAinf; %steady state minimum value of "a" at the end of each off-phase
        
        for i=1:length(f)
            
            %on phase
            a1=@(t) A(i).*(1-exp(-t/taua))+minAinf(i).*exp(-t/taua);
            Ia1(i)=quadl(a1,0,d(i),tol);
            Fa1=@(t) F(A(i).*(1-exp(-t/taua))+minAinf(i).*exp(-t/taua));
            I1(i)=quadl(Fa1,0,d(i),tol);
            
            %off phase
            a2=@(t) maxAinf(i).*exp(-(t-d(i))/taua);
            Ia2(i)=quadl(a2,d(i),T(i),tol);
            Fa2=@(t) F(maxAinf(i).*exp(-(t-d(i))/taua));
            I2(i)=quadl(Fa2,d(i),T(i),tol);
            
        end
        
        %%% analytical computation of mean a %%%
        meana1A=(A.*(d-taua*(1-exp(-d/taua)))+minAinf*taua.*(1-exp(-d/taua)))./d;
        
        %%%% mean of b %%%% 
        meanb=10*(I1+I2)./T;
        
        %NOTE: the 10 is needed to reproduce exactly the figures in the
        %manuscript - this was a typo in the paper. However, the shapes of
        %the responses are unaffected; this simply changes the vertical
        %scaling of <b>. The b-equation corresponding to the figures would
        %be: taub*db/dt=10*F(a)-b
    end

end
