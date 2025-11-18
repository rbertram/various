% Fig4-5.m

% This computer code is from "Interpreting Frequency Responses to Dose-Conserved Pulsatile
% Input Signals in Simple Cell Signaling Motifs", by P. Fletcher, F. Clement, A. Vidal,
% J. Tabak, and R. Bertram, PLoS One, 9:e95613, 2014.

function fig4_5_b_model()
% This code reproduces the panels for figures 4 and 5 of
%
% Fletcher PA, Clement F, Vidal A, Tabak J, Bertram R (2014)
% Interpreting Frequency Responses to Dose-Conserved Pulsatile Input
% Signals in Simple Cell Signaling Motifs. PLoS ONE 9(4): e95613.
% doi:10.1371/journal.pone.0095613
%
% Figure 4. The optimal pulse shape depends on the shape of the
% nonlinearity (panels: figures 1-6)
%
% Figure 5. Duration compensated signals yield a flat frequency
% response in the b-model (panels: figures 7,8)


clear;

%Signals parameters
T0=100; d0=10; A0=1; f0=1/T0;
fmax=1/d0;


%% Sublinear

%Michaelis-Menten function
n=1; K=1; Kn=K^n;
F=@(x) x.^n./(x.^n+Kn); 

f=10.^linspace(-3,log10(fmax),100);
T=1./f; 

%amplitude
A=A0*f0./f; d=d0*ones(size(f));
T1=10^2; d1=d0; A1=A0*T1/T0; f1=1/T1;
T2=T1/2; d2=d0; A2=A0*T2/T0; f2=1/T2;

meanb=(F(A).*d+F(0).*(T-d))./T;

% plot <b> vs f
figure(1);
semilogx(A,meanb,'k','linewidth',2); hold on

T=T1; d=d1; A=A1;
meanb1=(F(A).*d+F(0).*(T-d))./T;
T=T2; d=d2; A=A2;
meanb2=(F(A).*d+F(0).*(T-d))./T;

semilogx(A1,meanb1,'o',A2,meanb2,'s','linewidth',2)
xlabel('A')
ylabel('<b>_{\infty}')
AX(1)=gca;
set(gca,'ylim',[0,0.1],'ytick',[0,0.05,0.1],'yticklabel',{'0','0.05','0.1'})
set(gca,'xlim',[10^(-1),10^1],'xtick',[0.1,1,10],'xticklabel',{'0.1','1','10'})
axis square

% plot F(A)
figure(2)

AA=0:0.01:2;
plot(AA,F(AA),'k','linewidth',2); hold on
plot(A1,F(A1),'o',A2,F(A2),'s','linewidth',2)
xlabel('A')
ylabel('F(A)')
AX(2)=gca;
set(gca,'ylim',[0,0.8],'ytick',0:0.4:0.8,'yticklabel',{'0','0.4','0.8'})
set(gca,'xlim',[0,1.5],'xtick',0:0.5:1.5,'xticklabel',{'0','0.5','1','1.5'})
axis square

%% Superlinear

%superlinear
n=2; a=1; %use sign of a for increasing/decreasing
F=@(x) a*x.^n;


f=10.^linspace(-3,log10(fmax),100);
T=1./f; 

%amplitude
A=A0*f0./f; d=d0*ones(size(f));
T1=10^2; d1=d0; A1=A0*T1/T0; f1=1/T1;
T2=T1/2; d2=d0; A2=A0*T2/T0; f2=1/T2;

meanb=(F(A).*d+F(0).*(T-d))./T;

% plot <b> vs f
figure(3)

semilogx(A,meanb,'k','linewidth',2); hold on

T=T1; d=d1; A=A1;
meanb1=(F(A).*d+F(0).*(T-d))./T;
T=T2; d=d2; A=A2;
meanb2=(F(A).*d+F(0).*(T-d))./T;

semilogx(A1,meanb1,'o',A2,meanb2,'s','linewidth',2)
xlabel('A')
ylabel('<b>_{\infty}')
AX(3)=gca;
set(gca,'ylim',[0,0.3],'ytick',0:0.1:0.3,'yticklabel',{'0','0.1','0.2','0.3'})
set(gca,'xlim',[10^(-1),10^1],'xtick',[0.1,1,10],'xticklabel',{'0.1','1','10'})
axis square

% plot F(A) %%%%%%%%%%%%%%
figure(4)

AA=0:0.01:2;
plot(AA,F(AA),'k','linewidth',2); hold on
plot(A1,F(A1),'o',A2,F(A2),'s','linewidth',2)
xlabel('A')
ylabel('F(A)')
AX(4)=gca;
set(gca,'ylim',[0,1.5],'ytick',0:0.5:1.5,'yticklabel',{'0','0.5','1','1.5'})
set(gca,'xlim',[0,1.5],'xtick',0:0.5:1.5,'xticklabel',{'0','0.5','1','1.5'})
axis square


%% Sigmoid

%amplitude
f=10.^linspace(-3,log10(fmax),100);
T=1./f; 
A=A0*f0./f; d=d0*ones(size(f));

%hill function
n=2; K=A0;
F=@(x) x.^n./(x.^n+K^n); 

meanb3=(F(A).*d+F(0).*(T-d))./T;

%some points:
Mult=2;
A1=Mult*K; T1=Mult*T0; d1=d0; f1=1/T1;
A2=K; T2=T0;  d2=d0; f2=1/T2;
A3=K/Mult; T3=T0/Mult; d3=d0;  f3=1/T3;

mean1=(F(A1).*d1+F(0).*(T1-d1))./T1;
mean2=(F(A2).*d2+F(0).*(T2-d2))./T2;
mean3=(F(A3).*d3+F(0).*(T3-d3))./T3;

figure(5)

semilogx(A,meanb3,'k','linewidth',2);hold on
plot(A1,mean1,'bo','linewidth',2)
plot(A2,mean2,'gs','linewidth',2,'color',[0,0.5,0]);
plot(A3,mean3,'rd','linewidth',2)

xlabel('A');
ylabel('<b>_{\infty}');
AX(5)=gca;
set(gca,'ylim',[0,0.06],'ytick',0:0.02:0.06,'yticklabel',{'0','0.02','0.04','0.06'})
set(gca,'xlim',[10^(-1),10^1],'xtick',[0.1,1,10],'xticklabel',{'0.1','1','10'})
axis square


% sigmoidal y-intercept
Fp=@(x) n*Kn.*x.^(n-1)./(Kn+x.^n).^2;
Fpp=@(x) n*Kn.*x.^(n-2).*((n-1)*Kn-(n+1)*x.^n)./(Kn+x.^n).^3;

figure(6)

A=0:0.01:3; maxA=max(A);
y1=F(A);
plot(A,y1,'k','linewidth',2);hold on

x0=K/Mult;
b0=F(x0)-Fp(x0)*x0;
y2=Fp(x0)*A+b0; 
plot(x0,F(x0),'rd','linewidth',2);
plot(A,y2,'r--','linewidth',1);
plot(0,b0,'rd','linewidth',2)
plot([0,maxA],[b0,b0],'r:','linewidth',1)

x0=K;
b0=F(x0)-Fp(x0)*x0;
y2=Fp(x0)*A+b0; 
plot(x0,F(x0),'gs','linewidth',2,'color',[0,0.5,0]);
plot(A,y2,'g--','linewidth',1,'color',[0,0.5,0]);
plot(0,b0,'gs','linewidth',2,'color',[0,0.5,0]);
plot([0,maxA],[b0,b0],'g:','linewidth',1,'color',[0,0.5,0]);

x0=Mult*K;
b0=F(x0)-Fp(x0)*x0;
y2=Fp(x0)*A+b0; 
plot(x0,F(x0),'bo','linewidth',2);
plot(A,y2,'b--','linewidth',1);
plot(0,b0,'bo','linewidth',2)
plot([0,maxA],[b0,b0],'b:','linewidth',1)

set(gca,'ylim',[-.15,1])
set(gca,'ytick',0:0.5:1)
xlabel('A');
ylabel('F(A)')

AX(6)=gca;
set(gca,'xlim',[0,3])
set(gca,'ylim',[-0.2,1],'ytick',0:0.5:1,'yticklabel',{'0','0.5','1'})
axis square


%% Show that duration compensation doesn't work - Sigmoid

%amplitude
f=10.^linspace(-3,-1,100);
T=1./f; 
A=A0*ones(size(f)); d=d0*f0./f;

%hill function
n=2; K=1;
F=@(x) x.^n./(x.^n+K^n); 

meanb3=(F(A0).*d+F(0).*(T-d))./T;

%one point
T1=200; d1=d0*T1/T0; f1=1/T1; A1=A0;
T2=50; d2=d0*T2/T0; f2=1/T2; A2=A0;

mean1=(F(A0).*d1+F(0).*(T2-d1))./T1;
mean2=(F(A0).*d2+F(0).*(T2-d2))./T2;

AA=0:0.01:4;
FA=F(AA);

% subplot(1,2,1)
figure(7)

semilogx(f,meanb3,'k','linewidth',2);hold on
plot(f1,mean1,'o','linewidth',2)
plot(f2,mean2,'s','color',[0,0.5,0],'linewidth',2)
xlabel('f');
ylabel('<b>_{\infty}');
AX(7)=gca;
set(gca,'ylim',[0,0.1],'ytick',0:0.05:0.1,'yticklabel',{'0','0.05','0.1'})
set(gca,'xlim',[1e-3,1e-1],'xtick',[0.001,0.01,0.1],'xticklabel',{'0.001','0.01','0.1'})
axis square

% subplot(1,2,2)
figure(8)

plot(AA,FA,'k','linewidth',2);hold on
plot([K;K],[F(K);F(K)],'bo','linewidth',2)
plot([K;K],[F(K);F(K)],'s','color',[0,0.5,0],'linewidth',2)
xlabel('A');
ylabel('F(A)');
AX(8)=gca;
set(gca,'xlim',[0,2])
set(gca,'ylim',[0,1],'ytick',0:0.5:1,'yticklabel',{'0','0.5','1'})
axis square

%% format and save as eps
% 
% set(AX,'fontsize',18) %axis ticks
% set(findall(AX,'type','text'),'fontsize',18) %axis labels, etc.
% 
% FIG=1:8;
% set(FIG,'units','inches','position',[0.25,0.25,4,4]);
% % set(FIG,'paperunits','inches','paperposition',[0.25,0.25,4,4]);
% % set(FIG,'units','inches','outerposition',[0.2,0.2,4,4]);
% set(FIG,'color','w','paperpositionmode','auto');
% 
% POS=get(AX(1),'position'); %get axis dimensions of figure 1
% set(AX,'position',POS) %set all axes to same dimensions
% 
% fnames={'sublin_bA','sublin_FA','superlin_bA','superlin_FA',...
%     'sigmoid_bA','sigmoid_FA','sigmoid_dur_bf','sigmoid_dur_FA'};
% 
% for i=1:8
% %     print(i,'-depsc2','-painters',fnames{i});  %less white space
%     print(i,'-depsc2','-painters','-loose',fnames{i}); %more white space (full 4x4 in)
% end
%
% fix_lines('sigmoid_FA.eps')
