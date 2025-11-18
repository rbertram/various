% Fig2-3.m

% This computer code is from "Interpreting Frequency Responses to Dose-Conserved Pulsatile
% Input Signals in Simple Cell Signaling Motifs", by P. Fletcher, F. Clement, A. Vidal,
% J. Tabak, and R. Bertram, PLoS One, 9:e95613, 2014.

function fig2_3_minimalReceptor()
% This code reproduces the panels for figures 2 and 3 of
%
% Fletcher PA, Clement F, Vidal A, Tabak J, Bertram R (2014)
% Interpreting Frequency Responses to Dose-Conserved Pulsatile Input
% Signals in Simple Cell Signaling Motifs. PLoS ONE 9(4): e95613.
% doi:10.1371/journal.pone.0095613
%
% Figure 2. Responses of the ligand-receptor binding model to dose
% conserved input pulses. (panels: figures 1-3)
%
% Figure 3. The frequency response differs for systems with
% different steady state nonlinearity (panels: figures 4-7)

%parameters
p.gam0=0; 
p.k1p=.5; p.k1m=.25; K=p.k1m/p.k1p; %ligand binding
p.k2p=.1; p.k2m=.1; K2=p.k2m/p.k2p; %receptor dimerization

%The reference input signal parameters
T0=100; A0=1; d0=10; f0=1/T0; beta0=d0/T0;

%%  TIMECOURSES - different modes of compensation
nT0=2; fMult=4;
T1=T0; d1=d0; A1=A0; F1=1/T1;
T2=T1/fMult; d2=d0; A2=A0*T2/T0; F2=1/T2; %amplitude-frequency
T3=T1/fMult; d3=d0*T3/T0; A3=A0; F3=1/T3; %duration-frequency
T4=T1; d4=d0; A4=A0; F4=1/T4; %amplitude-duration

%time points to plot
dt=0.05;
tf=nT0*T0;
t1=0:dt:tf;

% Low frequency point %%
T=T0; d=d0; A=A0;
nT=nT0; %number of periods to plot

[L1,RL1,meanRL1]=ligandBindingTimecourse();

figure(1);clf

line([0,nT0*T1],[meanRL1,meanRL1],'color',[0.5,0,0.5],'linewidth',2)
line(t1,L1,'color','k','linewidth',1);
line(t1,RL1,'color',[0.5,0,0.5],'linewidth',2);

set(gca,'xlim',[0,nT0*T1-1.5])
box off
xlabel('t')
ylabel('RL')
axis square
ax1Y=get(gca,'ylim');

% High frequency point - Amplitude %%

T=T2; d=d2; A=A2;
nT=nT0*fMult;

[L2,RL2,meanRL2]=ligandBindingTimecourse();

figure(2);clf

line([0,nT0*fMult*T2],[meanRL2,meanRL2],'color',[0,0.5,0],'linewidth',2)
line(t1,L2,'color','r','linewidth',1);
line(t1,RL2,'color',[0,0.5,0],'linewidth',2);

axis square
box off
xlabel('t')
ylabel('RL')
set(gca,'xlim',[0,nT0*fMult*T2-1.5],'ylim',ax1Y);


% High frequency point - Duration %%

T=T3; d=d3; A=A3;
nT=nT0*fMult;

[L3,RL3,meanRL3]=ligandBindingTimecourse();

figure(3);clf

line([0,nT0*fMult*T3],[meanRL3,meanRL3],'color',[0,0.5,0],'linewidth',2)
line(t1,L3,'color','b','linewidth',1);
line(t1,RL3,'color',[0,0.5,0],'linewidth',2);

axis square
box off
xlabel('t')
ylabel('RL')
set(gca,'xlim',[0,nT0*fMult*T3-1.5],'ylim',ax1Y);



%% FREQUENCY RESPONSE 

% Frequency response for ligand binding is computed by integrating the
% analyical solution for the response to pulses.

tol=1e-6; %tolerance for quadl integration

%amplitude correction
f=10.^(-3:0.005:-1);
T=1./f;  
A=A0*T./T0; 
d=d0*ones(size(f));


for i=1:length(f)
    binf=A(i)/(A(i)+K);
    delon=p.k1p*A(i)+p.k1m;
    alph=binf*(1-exp(-delon*d(i)));
    bmaxinfAMP(i)=alph./(1-exp(-delon*d(i))*exp(-p.k1m*(T(i)-d(i))));
    bmininfAMP(i)=bmaxinfAMP(i)*exp(-p.k1m*(T(i)-d(i)));

    RL1=@(t) binf.*(1-exp(-delon.*t))+bmininfAMP(i).*exp(-delon.*t);
    I1(i)=quadl(RL1,0,d(i),tol);
    
    RL2=@(t) bmaxinfAMP(i).*exp(-p.k1m.*t);
    I2(i)=quadl(RL2,0,T(i)-d(i),tol);
end

f1=f;

meanRLamp=f.*(I1+I2);


%duration correction
f=10.^(-3:0.005:0); T=1./f;
A=A0*ones(size(f)); 
d=d0*T./T0;


for i=1:length(f)
    binf=A(i)/(A(i)+K);
    delon=p.k1p*A(i)+p.k1m;
    alph=binf*(1-exp(-delon*d(i)));
    bmaxinfDUR(i)=alph./(1-exp(-delon*d(i))*exp(-p.k1m*(T(i)-d(i))));
    bmininfDUR(i)=bmaxinfDUR(i)*exp(-p.k1m*(T(i)-d(i)));

    RL1=@(t) binf.*(1-exp(-delon.*t))+bmininfDUR(i).*exp(-delon.*t);
    I1(i)=quadl(RL1,0,d(i),tol);
    
    RL2=@(t) bmaxinfDUR(i).*exp(-p.k1m.*t);
    I2(i)=quadl(RL2,0,T(i)-d(i),tol);
end

f2=f;
meanRLdur=f.*(I1+I2);


figure(4);clf

semilogx(f1,meanRLamp,'r:','linewidth',2)
hold on
semilogx(f2,meanRLdur,'b--','linewidth',2)

ix1=find(f1>=F1,1,'first');
ix2=find(f1>=F2,1,'first');
ix3=find(f2>=F3,1,'first');

line(F1,meanRLamp(ix1),'color',[0.5,0,0.5],'marker','o','linewidth',2)
line(F2,meanRLamp(ix2),'color',[0,0.5,0],'marker','s','linewidth',2)
line(F3,meanRLdur(ix3),'color',[0,0.5,0],'marker','d','linewidth',2)

ylabel('<RL>_{\infty}');xlabel('f')
axis tight
axis square
set(gca,'ylim',[0,0.2],'xtick',[0.001,0.01,0.1,1],'xticklabel',{'0.001','0.01','0.1','1'})



%% Steady state nonlinearity, ligand binding
x=0:0.01:2;
RLss=x./(x+K);

figure(5);clf
plot(x,RLss,'k')
xlabel('L');ylabel('RL');
axis([0,2,0,1])
axis square




%% Dimerization

% Amplitude correction
df=0.25; %this sets the grid in frequency for dimerization frequency response

%Amplitude Compensation
f=10.^(-3:df:-1);
A=A0*f0./f; d=d0*ones(size(f));
T=1./f;


%Use numerical method to compute timecourse. First run the simulation for
%some time to remove initial transients, then average the solution values
%over a 5 periods

ttrans=ceil(1000./T).*T;
tf=5*T;
nF=length(f);

tic
for j=1:nF
    p.T=T(j); p.gam1=A(j); p.t1=d(j); disp(['f=' num2str(f(j)) ' A=' num2str(p.gam1) ' d=' num2str(p.t1)]);
    odeopt=odeset('RelTol',1e-6,'AbsTol',1e-5,'maxstep',p.t1/10);

    %remove transient
    yo=[0;0;0];
    [t,Y]=ode15s(@(t,Y) minRdimer(t,Y,p),[0,ttrans(j)],yo,odeopt);
    yo=[Y(end,1:2)';0];

    [t,Y]=ode15s(@(t,Y) minRdimer(t,Y,p),[0,tf(j)],yo,odeopt);
    DINT(j)=Y(end,3);
    meanD1(j)=DINT(j)/tf(j);
end
toc

f1=f;

% Duration Compensation
f=10.^(-3:df:0);
A=A0*ones(size(f)); d=d0*f0./f;
T=1./f;

ttrans=ceil(1000./T).*T;
tf=5*T;
nF=length(f);

%solve
tic
for j=1:nF
    p.T=T(j); p.gam1=A(j); p.t1=d(j); disp(['f=' num2str(f(j)) ' A=' num2str(p.gam1) ' d=' num2str(p.t1)]);
    odeopt=odeset('RelTol',1e-6,'AbsTol',1e-5,'maxstep',p.t1/10);

    %remove transient
    yo=[0;0;0];
    [t,Y]=ode15s(@(t,Y) minRdimer(t,Y,p),[0,ttrans(j)],yo,odeopt);
    yo=[Y(end,1:2)';0];

    [t,Y]=ode15s(@(t,Y) minRdimer(t,Y,p),[0,tf(j)],yo,odeopt);
    DINT(j)=Y(end,3);
    meanD2(j)=DINT(j)/tf(j);
end
f2=f;
toc

figure(6);
semilogx(f1',meanD1','r:','linewidth',2);hold on
semilogx(f2',meanD2','b--','linewidth',2)
ylabel('<D>_{\infty}');
xlabel('f')
set(gca,'xlim',[0.001,1],'ylim',[0,0.04],'xtick',[0.001,0.01,0.1,1],'xticklabel',{'0.001','0.01','0.1','1'})
axis square



%steady state nonlinearity, add dimerization
L=0:0.001:2;

%two branches..
x1=(-(1+K./L)+sqrt((1+K./L).^2+8/K2))./(4/K2);
x2=(-(1+K./L)-sqrt((1+K./L).^2+8/K2))./(4/K2);
y1=x1.^2/K2;
y2=x2.^2/K2;

figure(7)
plot(L,y1,'k')
xlabel('L'); ylabel('D')
axis([0,2,0,0.25])
axis square




%%%%% HELPER FUNCTIONS FOR COMPUTING TIMECOURSES AND MEANS

    function [L, RL, meanRL]=ligandBindingTimecourse()
        %t1, A, d, T, p.k1p, p.k1m are all from caller's workspace
        
        %input pulses
        L=zeros(1,length(t1));
        for i=1:length(t1)
            tmp=mod(t1(i),T); %phase of current time in [0,T) - L is on until dur
            L(i)=A*(tmp<=d);
        end

        %compute the exact steady state timecourse, one phase at a time
        ta=0:dt:d;
        tb=dt:dt:T-d;

        binf=A/(A+K);
        delon=p.k1p*A+p.k1m;
        alph=binf*(1-exp(-delon*d));
        maxXinf=alph/(1-exp(-delon*d)*exp(-p.k1m*(T-d)));
        minXinf=maxXinf*exp(-p.k1m*(T-d));

        RL=binf.*(1-exp(-delon.*ta))+minXinf.*exp(-delon.*ta);
        RL=[RL, maxXinf.*exp(-p.k1m.*tb)];

        ta=dt:dt:d;
        for i=2:nT
            RL=[RL, binf.*(1-exp(-delon.*ta))+minXinf.*exp(-delon.*ta)];
            RL=[RL, maxXinf.*exp(-p.k1m.*tb)];
        end

        meanRL=1/tf*trapz(t1,RL);
    end


end


function dY = minRdimer(t,Y,p)

    RL=Y(1);
    D=Y(2);

    %Generate the L value (square wave)
    tmp=mod(t,p.T); %phase of current time in [0,T) - L is on until p.dur
    L=p.gam1*(tmp<=p.t1)+p.gam0;
 
%     %Generate the L value (exponential drop)
%     tmp=mod(t,p.T); %phase of current time in [0,T) - L is on until p.dur
%     L=p.gam1*exp(-tmp)+p.gam0;
    
    r=1-RL-2*D;
    
    dRL = p.k1p*L*r - p.k1m*RL + 2*p.k2m*D -2*p.k2p*RL^2; 
    dD = p.k2p*RL^2 - p.k2m*D;
    da = D;
    
    dY=[dRL;dD;da;];
end
