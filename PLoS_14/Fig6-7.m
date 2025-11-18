% Fig6-7.m

% This computer code is from "Interpreting Frequency Responses to Dose-Conserved Pulsatile
% Input Signals in Simple Cell Signaling Motifs", by P. Fletcher, F. Clement, A. Vidal,
% J. Tabak, and R. Bertram, PLoS One, 9:e95613, 2014.

function fig6_7_ab_model()
% This code reproduces the panels for figures 6 and 7 of
%
% Fletcher PA, Clement F, Vidal A, Tabak J, Bertram R (2014)
% Interpreting Frequency Responses to Dose-Conserved Pulsatile Input
% Signals in Simple Cell Signaling Motifs. PLoS ONE 9(4): e95613.
% doi:10.1371/journal.pone.0095613
%
% Figure 6. Adding a linear component to the b-model results in frequency
% sensitivity to duration compensated signals (panels: figures 1-3)
%
% Figure 7. The ab-model displays the responses predicted by the nonlinear
% production rate function (panels: 4-6)

% NOTE: There is a minor typo in the paper regarding the b-equation. 
% The b-equation corresponding to the figures would
% be: taub*db/dt=10*F(a)-b.  
% The shapes of the responses are unaffected; this simply changes 
% the vertical scaling of <b>.  See line 212 of this code.

taua=10;

%Signals
T0=100; A0=3; beta0=0.1; f0=1/T0; d0=beta0*T0;

%numerical integral tolerance
tol=1e-6;

%% Figure 6 panels: How duration works in the ab-model

%duration correction
f=10.^(-4:0.01:-1); T=1./f;
A=A0*ones(size(f));
d=d0*f0./f;

%Choose two points to use for plotting timecourses of "a"
T1=200; d1=d0*T1/T0; f1=1/T1;
T2=50; d2=d0*T2/T0; f2=1/T2;

%analytical mean a and mean during on and off phases:
alph=A.*(1-exp(-d/taua));
maxAinf=alph./(1-exp(-T/taua));
minAinf=exp(-(T-d)/taua).*maxAinf;
ampInf=maxAinf-minAinf;


%calculate mean a ON/Off etc. for the two points.
alph=A0.*(1-exp(-[d1,d2]/taua));
maxAinfp=alph./(1-exp(-[T1,T2]/taua));
minAinfp=exp(-([T1,T2]-[d1,d2])/taua).*maxAinfp;

meanap=A0.*[d1,d2]./[T1,T2];
meanaONp=(A0.*([d1,d2]-taua*(1-exp(-[d1,d2]/taua)))+minAinfp*taua.*(1-exp(-[d1,d2]/taua)))./[d1,d2];
meanaOFFp=maxAinfp*taua.*(1-exp(-([T1,T2]-[d1,d2])/taua))./([T1,T2]-[d1,d2]);


% Some timecourses showing duration->amplitude, using above points
figure(1);clf

%T1
%first part: on phase
t1=linspace(0,d1,100);
a1=A0.*(1-exp(-t1/taua))+minAinfp(1).*exp(-t1/taua);
%second part: off phase
t2=linspace(d1,T1,100);
a2=maxAinfp(1).*exp(-(t2-d1)/taua);
    
plot([t1,t2],[a1,a2],'k','linewidth',2);hold on

%T2
%first part: on phase
t1=linspace(0,d2,100);
a1=A0.*(1-exp(-t1/taua))+minAinfp(2).*exp(-t1/taua);
%second part: off phase
t2=linspace(d2,T2,100);
a2=maxAinfp(2).*exp(-(t2-d2)/taua);
    
plot([t1,t2],[a1,a2],'k','linewidth',2);

set(gca,'xlim',[0,50])
xlabel('t');ylabel('a(t)')
axis square

%means
plot([0,d1],[meanaONp(1),meanaONp(1)],'-o','linewidth',2)
plot([0,d2],[meanaONp(2),meanaONp(2)],'-s','color',[0,0.5,0],'linewidth',2)


%% Now show results for a Hill function
n=2; K=1; Kn=K^n;
F=@(x) x.^n./(x.^n+Kn);

%F(a)
figure(2);clf
xx=0:0.01:2.5; plot(xx,F(xx),'k','linewidth',2);hold on
plot(meanaONp(1),F(meanaONp(1)),'o','linewidth',2);
plot(meanaONp(2),F(meanaONp(2)),'s','color',[0,0.5,0],'linewidth',2);
xlabel('a');ylabel('F(a)'); axis square
    
%<b>
figure(3);clf
%duration correction
A=A0*ones(size(f));
d=d0*f0./f;
meanb=computeMeans();

semilogx(f,meanb,'k','linewidth',2);hold on
xlabel('f');ylabel('<b>_{\infty}')
set(gca,'xlim',[0.0001,0.1]);
set(gca,'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
axis square

ix1=find(f>=f1,1,'first');
ix2=find(f>=f2,1,'first');
plot([f1,f1],[meanb(ix1),meanb(ix1)],'-o','linewidth',2)
plot([f2,f2],[meanb(ix2),meanb(ix2)],'-s','color',[0,0.5,0],'linewidth',2)


%% Figure 7: Frequency responses in the a-b model

% SUBLINEAR: Michaelis-Menten function
n=1; K=1; Kn=K^n;
F=@(x) x.^n./(x.^n+Kn);

%amplitude correction
A=A0*f0./f;
d=d0*ones(size(f));
meanb=computeMeans();

figure(4);clf
semilogx(f,meanb,'r:','linewidth',2);hold on

%duration correction
A=A0*ones(size(f));
d=d0*f0./f;
meanb=computeMeans();

semilogx(f,meanb,'b--','linewidth',2)
xlabel('f');ylabel('<b>_{\infty}')
set(gca,'xlim',[0.0001,0.1],'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
axis square

%SUPERLINEAR: power function
n=2; a=1; %use sign of a for increasing/decreasing
F=@(x) a*x.^n;

%amplitude correction
A=A0*f0./f;
d=d0*ones(size(f));
meanb=computeMeans();

figure(5);clf
loglog(f,meanb,'r:','linewidth',2);hold on

%duration correction
A=A0*ones(size(f));
d=d0*f0./f;
meanb=computeMeans();

loglog(f,meanb,'b--','linewidth',2)
xlabel('f');ylabel('<b>_{\infty}')
set(gca,'xlim',[0.0001,0.1],'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
set(gca,'ytick',[0.1,1,10,100],'yticklabel',{'0.1','1','10','100'})
axis square


%SIGMOIDAL: hill function
n=2; K=1; Kn=K^n;
F=@(x) x.^n./(x.^n+Kn);

%amplitude correction
A=A0*f0./f;
d=d0*ones(size(f));
meanb=computeMeans();

figure(6);clf
semilogx(f,meanb,'r:','linewidth',2);hold on

%duration correction
A=A0*ones(size(f));
d=d0*f0./f;
meanb=computeMeans();

semilogx(f,meanb,'b--','linewidth',2)
xlabel('f');ylabel('<b>_{\infty}')
set(gca,'xlim',[0.0001,0.1]);
set(gca,'xtick',[0.0001,0.001,0.01,0.1],'xticklabel',{'0.0001','0.001','0.01','0.1'})
axis square


%%% NESTED HELPER FUNCTIONS
    function meanb=computeMeans()
        
        alph=A.*(1-exp(-d/taua));
        maxAinf=alph./(1-exp(-T/taua)); %steady state maximum value of "a" at the end of each on-phase
        minAinf=exp(-(T-d)/taua).*maxAinf; %steady state minimum value of "a" at the end of each off-phase
        
        for i=1:length(f)
            %compute <F(a)> for each of on and off phase
            
            %on phase 
            Fa1=@(t) F(A(i).*(1-exp(-t/taua))+minAinf(i).*exp(-t/taua));
            I1(i)=quadl(Fa1,0,d(i),tol);
            
            %off phase
            Fa2=@(t) F(maxAinf(i).*exp(-(t-d(i))/taua));
            I2(i)=quadl(Fa2,d(i),T(i),tol);
            
        end
        
        %%%% mean of b %%%%
        meanb=10*(I1+I2)./T;
        
        %NOTE: the 10 is needed to reproduce exactly the figures in the
        %manuscript - this was a typo in the paper. However, the shapes of
        %the responses are unaffected; this simply changes the vertical
        %scaling of <b>. The b-equation corresponding to the figures would
        %be: taub*db/dt=10*F(a)-b
        
    end

end
