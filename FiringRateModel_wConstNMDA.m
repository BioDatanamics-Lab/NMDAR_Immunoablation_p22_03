function [Max,Min,HW]=FiringRateModel_wConstNMDA(S,EIratio,noise)

% 2021-11-01

% Firing Rate Model: How does increased firing rate affect spatial
% representation
% Pablo Jercog
% 
% Based on FiringRateModels.m (in laburo/papiros/paper16_4/untitled folder/
% GraphWorkspace paper16_4/Firing Rate Models Oscillatory Inputs). See also
% Striatum Veronica/Striatum MSN Firing Rate
% Models/StriatumNetworkMSN_FiringRatemodel.m
%
% Input-output (I-O) functions include Sigmoid, Threshold-Linear,
% Threshold-Nonlinear and Smooth Threshold-Nonlinear

 %clearvars;
% close all;

KSE = 2;
        % 1:   Minimal model: 2 CA3-PYR, 1 CA1-INT, 1 CA1-PYR. All cells
        %      are 2D including an adaptation process (for generality). All
        %      cells are 2D including an adaptation process for generality)
        % 2:   Extended model: N CA3-PYR, M CA1-INT, 1 CA1-PYR. 
        % 101: 1D single cell (ThrNlin, I-O function)
        % 102: 2D single cell (adaptation, ThrNlin I-O function)
       

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightsalmon = [0.9 0.5 0.4];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];
mediumacquamarine = [0.4 0.8 0.6];

% Functions
heaviside=@(t) 0.5*(t == 0)+(t > 0);
Fsin=@(t,f,Ain) sin(2*pi*f*t/1000);
Fsinpos=@(t,f,Ain) (1+sin(2*pi*f*t/1000-pi/2))/2;
%Falpha=@(t,alpha,ti) heaviside(t-ti).*(t-ti).*exp(-alpha*(t-ti));
Falpha=@(t,taur,taud,ti) heaviside(t-ti).*(exp(-(t-ti)/taud)-exp(-(t-ti)/taur));
Fgauss=@(t,a,ti,sgma) a*exp(-(t-ti).^2/(2*sgma^2));
% Threshold Linear 
Thrshldlin=@(x) (x > 0).*x;
% Sigmoid I-O function
xth = 0;
xslp = 1;
Sigexp=@(x) 1./(1+exp(-(x-xth)/xslp));
Invsigexp=@(x) xth+xslp*log(x./(1-x));
Sigexp1=@(x) 1./(1+exp(-(x-xth)/xslp)).^2*exp(-(x-xth)/xslp)/xslp;
% Threshold Nonlinear I-O function
Fthrnlin=@(x) 0*heaviside(-x)+heaviside(x).*heaviside(1-x).*x.^2+heaviside(x-1).*2.*sqrt(x-3/4);
Fthrnlin1=@(x) 0*heaviside(-x)+heaviside(x).*heaviside(1-x).*2.*x+heaviside(x-1).*1./sqrt(x-3/4);
InvFthrnlin=@(x) heaviside(x).*heaviside(1-x).*sqrt(x)+heaviside(x-1).*(x.^2+3)/4;
% Adaptation functions
Ainf=@(x,xth,xslp) 1./(1+exp(-(x-xth)/xslp));
InvAinf=@(x,xth,xslp) xth+xslp*log(x./(1-x));
Ainf1=@(x,xth,xslp) 1./(1+exp(-(x-xth)/xslp)).^2*exp(-(x-xth)/xslp)/xslp;
% Modified threshold Nonlinear I-O function 
% The direct definition of functions does not always work due to issues
% with the heaviside function. See specific (external) functions 
% SmoothThresholdNonlinear.m (Ftype=1, F IO function, Ftype=2, F', Ftype=3, Inv F)

if KSE == 1
    
    % Minimal model: 2 CA3-PYR, 1 CA1-INT, 1 CA1-PYR. All cells are 2D,
    % including an adaptation process (for generality).
    
    % Notation: CA3-PYR: Layer 1
    %           CA1-INT: Layer 2
    %           CA1-PYR: Layer 3
        
    % parameters 
    
        % Constant stimulation (inside the I-O function) 
        
    s1(1) = 0;
    s1(2) = 0;
    s2 = 0;
    s3 = 0;
    
        % Time constants
        
    taux1 = 20;
    taux2 = 10;
    taux3 = 20;
    tauz1 = 1;
    tauz2 = 1;
    tauz3 = 1;
    
        % Adaptation weight 
    
    g1(1) = 0;
    g1(2) = 0;
    g2 = 0;
    g3 = 0;
    
        % Stimulation (outside the I-O function)
    
    I1(1) = 0;
    I1(2) = 0;
    I2 = 0;
    I3 = 0;
    
        % Noise variance
    
    D1(1) = 0;
    D1(2) = 0;
    D2 = 0;
    D3 = 0;
    
        % Adaptation Ainf parameters
        
    ath = 0.2;
    aslp = 0.2;
    
        % Connectivity weights
        
        % Jba(n,m): connectivity from layer a to layer b
        %           m: cell number in layer a
        %           n: cell number in layer b
        
    J21(1,1) = 0.5;
    J21(1,2) = 0.5;
    J31(1,1) = 0.5;
    J31(1,2) = 0.5;
    J32(1,1) = 0.2;
          
    % Time definitions
    
    Tmax = 1000;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    % White noise

    eta = randn(1,length(t));
    
    % Place fields 
    
    ti(1) = 200;
    ti(2) = 300;
%     taur = 100;
%     taud = 101;
%     F(1,:) = Falpha(t,taur,taud,ti(1));
%     F(2,:) = Falpha(t,taur,taud,ti(2));

    a = 1;
    sgma = 50;
    F(1,:) = Fgauss(t,a,ti(1),sgma);
    F(2,:) = Fgauss(t,a,ti(2),sgma);
    F(1,:) = F(1,:)/max(F(1,:));
    F(2,:) = F(2,:)/max(F(2,:));
    
    % Numerical solution (modified Euler methods: Runge-Kutta 2)
    
        % Xa(n): Firing rate: layer a, cell n.
        % Za(n): Adaptation variable: layer a, cell n
        
    X1 = zeros(2,length(t));
    X2 = zeros(1,length(t));
    X3 = zeros(1,length(t));
    Z1 = zeros(2,length(t));
    Z2 = zeros(1,length(t));
    Z3 = zeros(1,length(t));
    
    for j=1:length(t)-1
        Arg = 0;
        k1x1(1) = (-X1(1,j)+SmoothThresholdNonlinear(Arg-g1(1)*Z1(1,j)+s1(1)+F(1,j),1)+I1(1))/taux1;
        Arg = 0;
        k1x1(2) = (-X1(2,j)+SmoothThresholdNonlinear(Arg-g1(2)*Z1(2,j)+s1(2)+F(2,j),1)+I1(2))/taux1;
        k1z1(1) = (Ainf(X1(1,j),ath,aslp)-Z1(1,j))/tauz1;
        k1z1(2) = (Ainf(X1(2,j),ath,aslp)-Z1(2,j))/tauz1;
        Arg = J21(1,1)*X1(1,j)+J21(1,2)*X1(2,j);
        k1x2(1) = (-X2(1,j)+SmoothThresholdNonlinear(Arg-g2(1)*Z2(1,j)+s2(1),1)+I2(1))/taux2;
        k1z2(1) = (Ainf(X2(1,j),ath,aslp)-Z2(1,j))/tauz2;
        Arg = J31(1,1)*X1(1,j)+J31(1,2)*X1(2,j)-J32(1,1)*X2(1,j);
        k1x3(1) = (-X3(1,j)+SmoothThresholdNonlinear(Arg-g3(1)*Z3(1,j)+s3(1),1)+I3(1))/taux3;    
        k1z3(1) = (Ainf(X3(1,j),ath,aslp)-Z3(1,j))/tauz3;       
        ax1(1) = X1(1)+k1x1(1)*dt;
        ax1(2) = X1(2)+k1x1(2)*dt;
        ax2(1) = X2(1)+k1x2(1)*dt;
        ax3(1) = X3(1)+k1x3(1)*dt;
        az1(1) = Z1(1)+k1z1(1)*dt;
        az1(2) = Z1(2)+k1z1(2)*dt;
        az2(1) = Z2(1)+k1z2(1)*dt;
        az3(1) = Z3(1)+k1z3(1)*dt;        
        Arg = 0;
        k2x1(1) = (-ax1(1)+SmoothThresholdNonlinear(Arg-g1(1)*az1(1)+s1(1)+F(1,j+1),1)+I1(1))/taux1;
        Arg = 0;
        k2x1(2) = (-ax1(2)+SmoothThresholdNonlinear(Arg-g1(2)*az1(2)+s1(2)+F(2,j+1),1)+I1(2))/taux1;
        k2z1(1) = (Ainf(ax1(1),ath,aslp)-az1(1))/tauz1;
        k2z1(2) = (Ainf(ax1(2),ath,aslp)-az1(2))/tauz1;
        Arg = J21(1,1)*ax1(1)+J21(1,2)*ax1(2);
        k2x2(1) = (-ax2(1)+SmoothThresholdNonlinear(Arg-g2(1)*az2(1)+s2(1),1)+I2(1))/taux2;
        k2z2(1) = (Ainf(ax2(1),ath,aslp)-az2(1))/tauz2;
        Arg = J31(1,1)*ax1(1)+J31(1,2)*ax1(2)-J32(1,1)*ax2(1);
        k2x3(1) = (-ax3(1)+SmoothThresholdNonlinear(Arg-g3(1)*az3(1)+s3(1),1)+I3(1))/taux3;    
        k2z3(1) = (Ainf(ax3(1),ath,aslp)-az3(1))/tauz3;        
        X1(1,j+1) = X1(1,j)+(k1x1(1)+k2x1(1))*dt/2;
        X1(2,j+1) = X1(2,j)+(k1x1(2)+k2x1(2))*dt/2;
        X2(1,j+1) = X2(1,j)+(k1x2(1)+k2x2(1))*dt/2;
        X3(1,j+1) = X3(1,j)+(k1x3(1)+k2x3(1))*dt/2;
        Z1(1,j+1) = Z1(1,j)+(k1z1(1)+k2z1(1))*dt/2;
        Z1(2,j+1) = Z1(2,j)+(k1z1(2)+k2z1(2))*dt/2;
        Z2(1,j+1) = Z2(1,j)+(k1z2(1)+k2z2(1))*dt/2;
        Z3(1,j+1) = Z3(1,j)+(k1z3(1)+k2z3(1))*dt/2;   
    end
    
    figure
    hold on
    plot(t,F(1,:),'-b','linewidth',2);
    plot(t,F(2,:),'-b','linewidth',2);
    set(gca,'fontsize',24);
    xlabel('t')
    ylabel('F');
    
    figure
    hold on
    plot(t,X1(1,:),'-b','linewidth',2);
    plot(t,X1(2,:),'-','Color',lightblueish,'linewidth',2);
    plot(t,X2(1,:),'-r','linewidth',2);
    plot(t,X3(1,:),'-g','linewidth',2);
     axis([0 Tmax 0 2]);
    set(gca,'fontsize',24);
    xlabel('t')
    ylabel('');
    legend('PYR1_{ca3}','PYR2_{ca3}','INT_{ca1}','PYR_{ca1}');
    
elseif KSE == 2
    
    % Extended model: M CA3-PYR, N CA1-INT, 1 CA1-PYR. All cells are 2D,
    % including an adaptation process (for generality).
    
    % Notation: CA3-PYR: Layer 1
    %           CA1-INT: Layer 2
    %           CA1-PYR: Layer 3
        
    % parameters 
    
    N = 21;% Number of CA3 pyramidal neurons
    M = 1;% Number of CA1 inhibitory neurons
    
        % Constant stimulation (inside the I-O function) 
            
    s1 = S*ones(1,N);%zeros(1,N);%1*ones(1,N);%zeros(1,N);
    s2 = S*ones(1,N);%zeros(1,M);
    s3 = S;
    
         % Constant stimulation (outside the I-O function) 
    %NMDA CONSTANT INPUT
    %NMDA=0.;%0.02;%0.05;%0.1;
    I1 = 0.0*ones(1,N);%zeros(1,N);
    I2 = 0.0*ones(1,N);%zeros(1,M);
    I3 = 0.0;%NMDA
        
        % Time constants
        
    taux1 = 20; %excitation time constant
    taux2 = 10; %inhibition time constant
    taux3 = 20; %excitation time constant 
    tauz1 = 1;%100 for example
    tauz2 = 1;
    tauz3 = 1;
    
        % Adaptation weight 
    
    g1 = 0.*ones(1,N);%zeros(1,N);
    g2 = 0.*ones(1,N);%zeros(1,M);
    g3 = 0.;
    
        % Noise variance
    
    D1 = noise*ones(1,N)/N;%zeros(1,N);
    D2 = noise*ones(1,M)/M;%zeros(1,M);
    D3 = noise;
    
    % Adaptation Ainf parameters
        
    ath = 0.2;
    aslp = 0.2;
        
    % NMDA parameters (NMDA "globally present") 
        
    Gnmdagl1 = zeros(1,N);%NMDA*ones(1,N);%zeros(1,N);
    Gnmdagl2 = zeros(1,M);%NMDA*ones(1,M);%zeros(1,M);
    Gnmdagl3 = 0;%NMDA;
    
         % Connectivity weights       
        % Jba(n,m): connectivity from layer a to layer b
        %           m: cell number in layer a
        %           n: cell number in layer b
       
    J21 = zeros(M,N);   
    J31 = zeros(1,N);
    J32 = zeros(1:M);
    INDINTENS=0.4;
    for k=1:N
        for l=1:M
            J21(l,k) = INDINTENS/N;
        end
    end
    for k=1:N
        J31(k) = EIratio*INDINTENS*Fgauss(k,1,N/2,5)/N;%gaussmf(k,[5 floor(N/2)])/N;%0.5; Fgauss=@(t,a,ti,sgma) a*exp(-(t-ti).^2/(2*sgma^2)); Fgauss(k,1,floor(N/2),5)
    end
    
    for l=1:M
        J32(l) = INDINTENS*4/M;
    end
    
     % Time definitions
    
    Tmax = 1000;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    % White noise

    %eta = randn(1,length(t));
   
    % LOCATION OF THE FIELDS
    
    ti = zeros(1,N);
    for k=1:N
       ti(k) = -(Tmax/(2*N))+k*(Tmax/N);
    end
    
    a = 1;
    sgma = 50;
    F = zeros(N,length(t));
    for k=1:N
        F(k,:) = Fgauss(t,a,ti(k),sgma);
        F(k,:) = F(k,:)/max(F(k,:));
    end
   
    % Numerical solution (modified Euler methods: Runge-Kutta 2)
    
        % Xa(n): Firing rate: layer a, cell n.
        % Za(n): Adaptation variable: layer a, cell n
        
    X1 = zeros(N,length(t));
    X2 = zeros(M,length(t));
    X3 = zeros(1,length(t));
    Z1 = zeros(N,length(t));
    Z2 = zeros(M,length(t));
    Z3 = zeros(1,length(t));
    k1x1 = zeros(1,N);
    k1z1 = zeros(1,N);
    k1x2 = zeros(1,M);
    k1z2 = zeros(1,M);
    ax1 = zeros(1,N);
    az1 = zeros(1,N);
    ax2 = zeros(1,M);
    az2 = zeros(1,M);
    k2x1 = zeros(1,N);
    k2z1 = zeros(1,N);
    k2x2 = zeros(1,M);
    k2z2 = zeros(1,M);
    
    for j=1:length(t)-1    
        eta1 = randn(1,N);
        eta2 = randn(1,M);
        eta3 = randn;
        for k=1:N
            Arg = 0;
            k1x1(k) = (-X1(k,j)+SmoothThresholdNonlinear(Arg-g1(k)*Z1(k,j)+s1(k)+Gnmdagl1(k)*X1(k,j)+F(k,j),1)+I1(k))/taux1;
            k1z1(k) = (Ainf(X1(k,j),ath,aslp)-Z1(k,j))/tauz1;
        end
        for l=1:M
            Arg = 0;
            for k=1:N
                Arg = Arg+J21(l,k)*X1(k,j);    
            end
            k1x2(l) = (-X2(l,j)+SmoothThresholdNonlinear(Arg-g2(l)*Z2(l,j)+s2(l)+Gnmdagl2(l)*X2(l,j),1)+I2(l))/taux2;
            k1z2(l) = (Ainf(X2(l,j),ath,aslp)-Z2(l,j))/tauz2;
        end
        Arg = 0;
        for k=1:N
            Arg = Arg+J31(1,k)*X1(k,j);
        end
        for l=1:M
            Arg = Arg-J32(1,l)*X2(l,j);
        end
        k1x3(1) = (-X3(1,j)+SmoothThresholdNonlinear(Arg-g3(1)*Z3(1,j)+s3(1)+Gnmdagl3*X3(1,j),1)+I3(1))/taux3;    
        k1z3(1) = (Ainf(X3(1,j),ath,aslp)-Z3(1,j))/tauz3;  
        for k=1:N
            ax1(k) = X1(k)+k1x1(k)*dt;
            ax1(k) = ax1(k)+sqrt(2*D1(k)*dt)*eta1(k);
            az1(k) = Z1(k)+k1z1(k)*dt;
        end
        for l=1:M
             ax2(l) = X2(l)+k1x2(l)*dt;
             ax2(l) = ax2(l)+sqrt(2*D2(l)*dt)*eta2(l);
             az2(l) = Z2(l)+k1z2(l)*dt;
        end
        ax3(1) = X3(1)+k1x3(1)*dt;
        ax3(1) = ax3(1)+sqrt(2*D3(1)*dt)*eta3(1);
        az3(1) = Z3(1)+k1z3(1)*dt;            
        for k=1:N
            Arg = 0;
            k2x1(k) = (-ax1(k)+SmoothThresholdNonlinear(Arg-g1(k)*az1(k)+s1(k)+Gnmdagl1(k)*ax1(k)+F(k,j),1)+I1(k))/taux1;
            k2z1(k) = (Ainf(ax1(k),ath,aslp)-az1(k))/tauz1;
        end
        for l=1:M
            Arg = 0;
            for k=1:N
                Arg = Arg+J21(l,k)*ax1(k);    
            end
            k2x2(l) = (-ax2(l)+SmoothThresholdNonlinear(Arg-g2(l)*az2(l)+s2(l)+Gnmdagl2(l)*ax2(l),1)+I2(l))/taux2;
            k2z2(l) = (Ainf(ax2(l),ath,aslp)-az2(l))/tauz2;
        end
        Arg = 0;
        for k=1:N
            Arg = Arg+J31(1,k)*ax1(k);
        end
        for l=1:M
            Arg = Arg-J32(1,l)*ax2(l);
        end
        k2x3(1) = (-ax3(1)+SmoothThresholdNonlinear(Arg-g3(1)*az3(1)+s3(1)+Gnmdagl3(1)*ax3(1),1)+I3(1))/taux3;    
        k2z3(1) = (Ainf(ax3(1),ath,aslp)-az3(1))/tauz3;  
        for k=1:N
            X1(k,j+1) = X1(k,j)+(k1x1(k)+k2x1(k))*dt/2;
            Z1(k,j+1) = Z1(k,j)+(k1z1(k)+k2z1(k))*dt/2;
            X1(k,j+1) = X1(k,j+1)+sqrt(2*D1(k)*dt)*eta1(k);
        end
        for l=1:M
            X2(l,j+1) = X2(l,j)+(k1x2(l)+k2x2(l))*dt/2;
            Z2(l,j+1) = Z2(l,j)+(k1z2(l)+k2z2(l))*dt/2;
            X2(l,j+1) = X2(l,j+1)+sqrt(2*D2(l)*dt)*eta2(l);
        end
        X3(1,j+1) = X3(1,j)+(k1x3(1)+k2x3(1))*dt/2;
        Z3(1,j+1) = Z3(1,j)+(k1z3(1)+k2z3(1))*dt/2; 
        X3(1,j+1) = X3(1,j+1)+sqrt(2*D3(1)*dt)*eta3(1);
    end
    
    figure(601);clf;
    
    subplot(2,2,1)
    hold on
    for k=1:N
        plot(t,F(k,:),'-b','linewidth',2);
    end
    set(gca,'fontsize',10);
    xlabel('t')
    ylabel('F');
    hold off
    
    subplot(2,2,2)
    hold on
    %for k=1:N
        plot(t,X3(1:end),'-g','linewidth',2);
    %end
    set(gca,'fontsize',10);
    xlabel('t')
    ylabel('F');
    %ylim([0,max(X3)*1.1])
    ylim([0,0.2])
    hold off
    
    subplot(2,2,3)
    hold on
        plot(t,X3(1:end),'-g','linewidth',2);
    for k=1:M
        plot(t,X2(k,1:end),'-r','linewidth',2);
    end
    set(gca,'fontsize',10);
    xlabel('t')
    ylabel('F');
    ylim([0,4])
    hold off
    
    subplot(2,2,4)
    hold on
    text(0,0.5,['Exc-In = ',num2str(Gnmdagl3*max(X3)+EIratio*INDINTENS*sum(Fgauss([1:1:N],1,N/2,5)/N))]);
    text(0,0.4,['Inh-In = ',num2str(Gnmdagl3*max(X2)+INDINTENS*4/M)]);
    text(0,0.3,['Exc/Inh (Max) = ',num2str((Gnmdagl3*max(X3)+EIratio*INDINTENS*sum(Fgauss([1:1:N],1,N/2,5)/N))/(Gnmdagl3*max(X2)+INDINTENS*4/M))]);
    text(0,0.2,['Exc/Inh (Min) = ',num2str((Gnmdagl3*min(X3)+EIratio*INDINTENS*sum(Fgauss([1:1:N],1,N/2,5)/N))/(Gnmdagl3*min(X2)+INDINTENS*4/M))]);
    axis off
    hold off
    
    Nframes=200;
    Window=(1/Nframes)*ones(1,Nframes);
    V=filtfilt(Window,1,X3);
    Max=max(V(1:end));
    Min=median(V(1:end));
    Mean=mean([Max,Min]);
    HW1=find(Mean<V(1:end),1,'first');
    HW2=find(Mean<V(1:end),1,'last');
    HW=((HW2-HW1)/length(V))*50;
    
elseif KSE == 101
    
     % 1D single cell (ThrNlin I-O function)
    
    % Parameters

    J = 1;
    taux = 1;
    s = 0;
    I = 0;
    D = 0;
    
    J = 0;

    % Time definitions
    
    Tmax = 20;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    % White noise

    eta = randn(1,length(t));
    
    % Numerical solution

    x = zeros(1,length(t));

    x(1) = 0;


    for j=1:length(t)-1
        kx1 = (-x(j)+SmoothThresholdNonlinear(J*x(j)-s,1)+I)/taux;
        ax = x(j)+kx1*dt;
        ax = ax+sqrt(2*D*dt)*eta(j);
        kx2 = (-ax+SmoothThresholdNonlinear(J*x(j)-s,1)+I)/taux;
        x(j+1) = x(j)+(kx1+kx2)*dt/2;    
        x(j+1) = x(j+1)+sqrt(2*D*dt)*eta(j);
    end
    
    figure
    hold on
    plot(t,x,'-b','linewidth',2);
    axis([0 Tmax 0 1]);
    set(gca,'fontsize',24);
    xlabel('t');
    ylabel('x');
    title('activity vs. time');
    
    xvec = -5:0.025:5;
    
    figure
    hold on
    plot(xvec,SmoothThresholdNonlinear(J*xvec-s,1)+I,'-b','linewidth',2);
    plot(xvec,xvec,'-r','linewidth',2);
    axis([-1 4 -0.2 4]);
    xlabel('x');
    ylabel('F(x)');
    set(gca,'fontsize',24);
    
    figure
    hold on
    plot(xvec,zeros(1,length(xvec)),'--','Color',[.7 .7 .7],'linewidth',1);
    plot(xvec,-xvec+SmoothThresholdNonlinear(J*xvec-s,1)+I,'-b','linewidth',2);
    axis([0 5 -1 1]);
    xlabel('x');
    ylabel('N_x(x)');
    set(gca,'fontsize',24);
    title('Phase-space diagram');
    
elseif KSE == 102 
    
    % 2D single cell (adaptaion, ThrNlin I-O function)
     
    % Parameters
    
    J = -0.1;
    tau = 1;
    s = 0.75;
    I = 0;
    g = 1;
    ath = 0.2;
    aslp = 0.2;
    
    taux = 1;
    tauz = 10;
    
    I = 0;
    D = 0;

    % Time definitions
    
    Tmax = 1000;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    % White noise

    eta = randn(1,length(t));
    
    % Numerical solution

    x = zeros(1,length(t));
    z = zeros(1,length(t));
    
    for j=1:length(t)-1
        k1x = (-x(j)+SmoothThresholdNonlinear(J*x(j)-g*z(j)+s,1))/taux;
        k1z = (Ainf(x(j),ath,aslp)-z(j))/tauz;
        ax = x(j)+k1x*dt;
        ax = ax+sqrt(2*D*dt)*eta(j);
        az = z(j)+k1z*dt;
        k2x = (-ax+SmoothThresholdNonlinear(J*ax-g*az+s,1))/taux;
        k2z = (Ainf(ax,ath,aslp)-az)/tauz;
        x(j+1) = x(j)+(k1x+k2x)*dt/2;
        x(j+1) = x(j+1)+sqrt(2*D*dt)*eta(j);
        z(j+1) = z(j)+(k1z+k2z)*dt/2;        
    end
    
    figure
    hold on
    plot(t,x,'-b','linewidth',2);
    plot(t,z,'-r','linewidth',2);
    set(gca,'fontsize',24);
    axis([0 Tmax 0 5]);
    xlabel('t');
    ylabel('');
    legend('x_1','x_2');
    
    % Computation of the nullclines
    
    xvec = 0:0.01:5;
    Nlcx = (J*xvec-SmoothThresholdNonlinear(xvec,3)+s)/g;
    Nlcz = Ainf(xvec,ath,aslp);
    
    rconstraint = (s+J*xvec)/g;
    
    % Computattion of the fixed-point
    
    xfp = zeros(1);
    zfp = zeros(1);
    xstab = zeros(1,2);
    stab = zeros(1);
    
    F = Nlcx-Nlcz;
    
    cnt = 0;
    if F(1) > 0
        for j=2:length(xvec)
            if F(j)<0 && cnt == 0
                cnt = cnt+1;
                xfp(cnt) = xvec(j);
                zfp(cnt) = Nlcz(j);
            elseif F(j)>0 && cnt == 1
                cnt = cnt+1;
                xfp(cnt) = xvec(j);
                zfp(cnt) = Nlcz(j);
            elseif F(j)<0 && cnt == 2
                cnt = cnt+1;
                xfp(cnt) = xvec(j);
                zfp(cnt) = Nlcz(j);
            end
        end
    elseif F(1) < 0
        for j=2:length(xvec)
            if F(j)>0 && cnt == 0
                cnt = cnt+1;
                xfp(cnt) = xvec(j);
                zfp(cnt) = Nlcz(j);
            elseif F(j)<0 && cnt == 1
                cnt = cnt+1;
                xfp(cnt) = xvec(j);
                zfp(cnt) = Nlcz(j);
            elseif F(j)>0 && cnt == 2
                cnt = cnt+1;
                xfp(cnt) = xvec(j);
                zfp(cnt) = Nlcz(j);
            end
        end
    end
    
     % Stability 
    
    a = (J*SmoothThresholdNonlinear(s+J*xfp-g*zfp,2)-1)/taux;
    b = -SmoothThresholdNonlinear(s+J*xfp-g*zfp,2)*g/taux;
    c = Ainf1(xfp,xth,xslp)/tauz;
    d = -1/tauz;
    
   
    
    r = roots([1 -(a+d) a*d-b*c]);
    if (a-d)^2+4*b*c<0
        mu = sqrt(-(4*b*c+(a-d)^2))/2;
    else
        mu = 0;
    end
    fnat = mu*1000/(2*pi);
    
   [real(r(1)) real(r(2)) fnat]
    
    hfig = figure;
    %figure
    hold on    
    plot(xvec,Nlcx,'-r','linewidth',2);
    plot(xvec,Nlcz,'-g','linewidth',2);
    plot(x,z,'-b','linewidth',2);
    set(gca,'fontsize',24);
    axis([0 5 0 3]);
    xlabel('x');
    ylabel('z');
    set(hfig, 'Position', [800 600 550 500]);
    plot(xvec,rconstraint,'--','Color',lightblueish,'linewidth',2);
    
    
end



        