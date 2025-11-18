% Fig10.m

% This computer code is from "Interpreting Frequency Responses to Dose-Conserved Pulsatile 
% Input Signals in Simple Cell Signaling Motifs", by P. Fletcher, F. Clement, A. Vidal, 
% J. Tabak, and R. Bertram, PLoS One, 9:e95613, 2014.

function fig10_phosphorylation()
% This code reproduces the panels for figure 10 of
%
% Fletcher PA, Clement F, Vidal A, Tabak J, Bertram R (2014)
% Interpreting Frequency Responses to Dose-Conserved Pulsatile Input
% Signals in Simple Cell Signaling Motifs. PLoS ONE 9(4): e95613.
% doi:10.1371/journal.pone.0095613
%
% Figure 10. Matching input signal characteristics to the
% sigmoidal nonlinearity and system kinetics in the Goldbeter-
% Koshland model

odeopt=odeset('RelTol',1e-8,'AbsTol',1e-6,'NonNegative',[1:4]);

%parameters
p.E2tot=50; p.Wtot=1000;
p.a1=50; p.d1=499; p.k1=1;
p.a2=50; p.d2=499; p.k2=1;
p.gam0=0;
lambda=1;

K1=(p.d1+ p.k1)/p.a1;
K2=(p.d2+ p.k2)/p.a2;

%===Signals===
T0=100; d0=10;
f=10.^(-3:0.1:-1); T=1./f;


%Amplitude compensation
d=d0*ones(size(f));

A0=100; A=A0*T./T0;
meanWs1=computeMean(p);

A0=350; A=A0*T./T0;
lambda=1;
meanWs2=computeMean(p);

lambda=0.2;
meanWs2b=computeMean(p);

lambda=5;
meanWs2c=computeMean(p);

lambda=1;
A0=600; A=A0*T./T0;
meanWs3=computeMean(p);

% save GK_A0_lambda_Amp

%plot amplitude compensation panels
%amp+A0
figure(1);clf
semilogx(f,meanWs1./p.Wtot,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
semilogx(f,meanWs2./p.Wtot,'-','color',[0,0.5,0],'linewidth',2)
semilogx(f,meanWs3./p.Wtot,'--','color',[0.5,0.5,0],'linewidth',2)
set(gca,'xtick',[0.001,0.01,0.1],'xticklabel',{'0.001','0.01','0.1'})
axis([0.001,0.1,0,0.9])
xlabel('f')
ylabel('<W^*>_{\infty}')
axis square

%amp350+lambda
figure(2);clf
semilogx(f,meanWs2b./p.Wtot,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
semilogx(f,meanWs2./p.Wtot,'-','color',[0,0.5,0],'linewidth',2)
semilogx(f,meanWs2c./p.Wtot,'--','color',[0.5,0.5,0],'linewidth',2)
set(gca,'xtick',[0.001,0.01,0.1],'xticklabel',{'0.001','0.01','0.1'})
axis([0.001,0.1,0,0.35])
xlabel('f');
ylabel('<W^*>_{\infty}')
axis square


%Duration compensation panels
d=d0*T./T0;

A0=100; A=A0*ones(size(f));
meanWs1=computeMean(p);

A0=350; A=A0*ones(size(f));
lambda=1;
meanWs2=computeMean(p);

lambda=0.2;
meanWs2b=computeMean(p);

lambda=5;
meanWs2c=computeMean(p);

lambda=1;
A0=600; A=A0*ones(size(f));
meanWs3=computeMean(p);

% save GK_A0_lambda_Dur

%plot duration compensation panels
%amp+A0
figure(3);clf
semilogx(f,meanWs1./p.Wtot,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
semilogx(f,meanWs2./p.Wtot,'-','color',[0,0.5,0],'linewidth',2)
semilogx(f,meanWs3./p.Wtot,'--','color',[0.5,0.5,0],'linewidth',2)
set(gca,'xtick',[0.001,0.01,0.1],'xticklabel',{'0.001','0.01','0.1'})
axis([0.001,0.1,0,0.9])
xlabel('f')
ylabel('<W^*>_{\infty}')
axis square

%amp350+lambda
figure(4);clf
semilogx(f,meanWs2b./p.Wtot,'-.','color',[0,0.5,0.5],'linewidth',2);hold on
semilogx(f,meanWs2./p.Wtot,'-','color',[0,0.5,0],'linewidth',2)
semilogx(f,meanWs2c./p.Wtot,'--','color',[0.5,0.5,0],'linewidth',2)
set(gca,'xtick',[0.001,0.01,0.1],'xticklabel',{'0.001','0.01','0.1'})
axis([0.001,0.1,0,0.35])
xlabel('f');
ylabel('<W^*>_{\infty}')
axis square

% NESTED FUNCTION

    function meanWs=computeMean(p)
        
        p.a1=p.a1*lambda; p.d1=p.d1*lambda; p.k1=p.k1*lambda; 
        p.a2=p.a2*lambda; p.d2=p.d2*lambda; p.k2=p.k2*lambda;
        
        nT=3;
        tf=nT.*T;
        
        tic
        for i=1:length(f)p.T=T(i); p.gam1=A(i); p.t1=d(i);
            disp(['f=' num2str(f(i)) ' A=' num2str(p.gam1) ' d=' num2str(p.t1)]);
            
            %initial conditions
            yo=[0;0;0;0];
            
            if d(i)==T(i)
                
                p.L=A(i);
                [t,Y]=ode15s(@(t,Y) GoldbeterKoshland(t,Y,p),[0,max(2*T(i),500)],yo,odeopt);
                yo=[Y(end,1);Y(end,2);Y(end,3);0];
                
                [t,Y]=ode15s(@(t,Y) GoldbeterKoshland(t,Y,p),[0,tf(i)],yo,odeopt);
                
            else
                
                %Transient: solve one phase at a time.
                t0=0;
                j=0;
                while t0<max(2*T(i),500)
                    p.L=A(i);
                    [t,Y]=ode15s(@(t,Y) GoldbeterKoshland(t,Y,p),[t0,d(i)+j*T(i)],yo,odeopt);
                    yo=[Y(end,1);Y(end,2);Y(end,3);0];
                    
                    p.L=0;
                    [t,Y]=ode15s(@(t,Y) GoldbeterKoshland(t,Y,p),j*T(i)+[d(i),T(i)],yo,odeopt);
                    yo=[Y(end,1);Y(end,2);Y(end,3);0];
                    
                    t0=t(end);
                    j=j+1;
                end
                
                
                %now solve the next few periods to measure the mean
                t0=0;
                j=0;
                while t0<tf(i)
                    p.L=A(i);
                    [t,Y]=ode15s(@(t,Y) GoldbeterKoshland(t,Y,p),[t0,d(i)+j*T(i)],yo,odeopt);
                    yo=Y(end,:);
                    
                    p.L=0;
                    [t,Y]=ode15s(@(t,Y) GoldbeterKoshland(t,Y,p),j*T(i)+[d(i),T(i)],yo,odeopt);
                    yo=Y(end,:);
                    
                    t0=t(end);
                    j=j+1;
                end
                
            end
            
            WINT(i)=Y(end,4);meanWs(i)=WINT(i)/tf(i);
            
            
        end
        toc
    end
end

function dY = GoldbeterKoshland(t,Y,p)

    Ws=Y(1);
    C1=Y(2);
    C2=Y(3);
    
    E1tot=p.L; %input changes total amount of active enzyme

    dWs= p.k1*C1 +p.d2*C2 -p.a2*Ws*(p.E2tot-C2);
    dC1= p.a1*(p.Wtot-Ws-C1-C2)*(E1tot-C1) -(p.d1+p.k1)*C1;
    dC2= p.a2*Ws*(p.E2tot-C2)-(p.d2+p.k2)*C2;

    dWsint = Ws;
    
    dY=[dWs;dC1;dC2;dWsint];
end
