alf=1/6; %Alpha parameter
bet=1/3000; %Beta parameter
k = 1/2; %k parameter
N=1000; %initial population size
pI=0.01; %fraction of population initially infected

S0=N.*(1-pI); %Initial Susceptible Population
E0=0; %Initial Exposed Population
I0=N.*pI; %Initial Infected Population
R0=0; %Initial Recovered Population

T=120; %time horizon
M=120; % number of time steps
Npaths=100; %number of sample paths
dt=T./M; % length of time increment
t=0:dt:T;

%Deterministic Model Initialization
St=zeros(1,length(t));
Et=zeros(1,length(t));
It=zeros(1,length(t));
Rt=zeros(1,length(t));

St(1)=S0;
Et(1)=E0;
It(1)=I0;
Rt(1)=R0;

%Stochastic Model Initialization
StNp=zeros(Npaths,length(t));
EtNp=zeros(Npaths,length(t));
ItNp=zeros(Npaths,length(t));
RtNp=zeros(Npaths,length(t));

StNp(:,1)=S0;
EtNp(:,1)=R0;
ItNp(:,1)=I0;
RtNp(:,1)=R0;

 for i=1:M
% Deterministic Model
    St(i+1)=St(i).*(1-bet.*It(i).*dt);
    Et(i+1)=bet.*St(i).*It(i).*dt+(1-k.*dt).*Et(i);
    It(i+1)=k.*Et(i).*dt+It(i).*(1-alf.*dt);
    Rt(i+1)=Rt(i)+alf.*It(i).*dt;

% Stochastic model
    dS=binornd(StNp(:,i),bet.*dt.*ItNp(:,i),Npaths,1);
    dE=binornd(EtNp(:,i),k.*dt,Npaths,1);
    dI=binornd(ItNp(:,i),alf.*dt,Npaths,1);
    
    StNp(:,i+1)=StNp(:,i)-dS;
    EtNp(:,i+1)=EtNp(:,i)+dS-dE;
    ItNp(:,i+1)=ItNp(:,i)+dE-dI;
    RtNp(:,i+1)=RtNp(:,i)+dI;
 end


% Solution to b
figure;
hold on
plot(t,St,'b','LineWidth',2)
plot(t,Et,'y','LineWidth',2)
plot(t,It,'g','LineWidth',2)
plot(t,Rt,'r','LineWidth',2)
xlabel('time')
ylabel('Population Group Size')
legend('Susceptibles','Exposed','Infecteds','Recovereds')
title('Deterministic Model')
hold off

figure;
hold on
plot(t,St,'b','LineWidth',2)
plot(t,Et,'y','LineWidth',2)
plot(t,It,'g','LineWidth',2)
plot(t,Rt,'r','LineWidth',2)

plot(t,StNp(1:Npaths,:),'b','LineWidth',1)
plot(t,EtNp(1:Npaths,:),'y','LineWidth',1)
plot(t,ItNp(1:Npaths,:),'g','LineWidth',1)
plot(t,RtNp(1:Npaths,:),'r','LineWidth',1)

plot(t,St,'k','LineWidth',2)
plot(t,Et,'k','LineWidth',2)
plot(t,It,'k','LineWidth',2)
plot(t,Rt,'k','LineWidth',2)

legend('Susceptibles','Exposed','Infecteds','Recovereds')

xlabel('time')
ylabel('Population Group Size')

title('Stochastic vs Deterministic')

hold off

% Solution to Part c
% Calculate mean and variance functions for S, I, R, and approx 95% CIs

AveSt=mean(StNp);
AveEt=mean(EtNp);
AveIt=mean(ItNp);
AveRt=mean(RtNp);

VSt=var(StNp);
VEt=var(EtNp);
VIt=var(ItNp);
VRt=var(RtNp);

UCLSt=AveSt+1.96.*sqrt(VSt./Npaths);
LCLSt=AveSt-1.96.*sqrt(VSt./Npaths);
UCLEt=AveEt+1.96.*sqrt(VEt./Npaths);
LCLEt=AveEt-1.96.*sqrt(VEt./Npaths);
UCLIt=AveIt+1.96.*sqrt(VIt./Npaths);
LCLIt=AveIt-1.96.*sqrt(VIt./Npaths);
UCLRt=AveRt+1.96.*sqrt(VRt./Npaths);
LCLRt=AveRt-1.96.*sqrt(VRt./Npaths);

figure;
hold on
plot(t,St,'b','LineWidth',2)
plot(t,Et,'y','LineWidth',2)
plot(t,It,'g','LineWidth',2)
plot(t,Rt,'r','LineWidth',2)

plot(t,AveSt,'k','LineWidth',1)
plot(t,AveEt,'k','LineWidth',1)
plot(t,AveIt,'k','LineWidth',1)
plot(t,AveRt,'k','LineWidth',1)

plot(t,UCLSt,':k','LineWidth',1)
plot(t,LCLSt,':k','LineWidth',1)
plot(t,UCLEt,':k','LineWidth',1)
plot(t,LCLEt,':k','LineWidth',1)
plot(t,UCLIt,':k','LineWidth',1)
plot(t,LCLIt,':k','LineWidth',1)
plot(t,UCLRt,':k','LineWidth',1)
plot(t,LCLRt,':k','LineWidth',1)
xlabel('time')
ylabel('Population Group Size')
legend('Susceptibles','Exposed','Infecteds','Recovereds')
title("Stochastic Model - Mean & CI")
hold off


%Solution to part d
[Imaximum T_Imaximum] = max(It)

%Solution to part e
figure;
histogram(ItNp(1:Npaths,T_Imaximum), 'Normalization','probability')
xlabel('Simulated Is')
ylabel('Probability')
ylim([0 0.6])
title("Simulated I's at the time of peak in Deterministic Model")

%Solution to part f
tao = find(It>0.9*Imaximum);
tao0 = tao(1)
tao1 = tao(length(tao))
length_stressed_period = tao1-tao0

%Solution to part g
taoNp = zeros(Npaths,3);
for i=1:Npaths
    ArrayTemp = find(ItNp(i,:)>0.9*Imaximum);
    if isempty(ArrayTemp)
        taoNp(i,1) = T;
        taoNp(i,2) = T;
        taoNp(i,3) = taoNp(i,2)-taoNp(i,1);
    else
        taoNp(i,1) = ArrayTemp(1);
        taoNp(i,2) = ArrayTemp(length(ArrayTemp));
        taoNp(i,3) = taoNp(i,2)-taoNp(i,1);
    end
end

tao0Np = mean(taoNp(1:Npaths,1))
tao1Np = mean(taoNp(1:Npaths,2))
taoStrPeriodNp = mean(taoNp(1:Npaths,3))

figure;
histogram(taoNp(1:Npaths,1), 'Normalization', 'probability')
xlabel('{\tau}_0')
ylabel('Probability')
ylim([0 1])
title('Histogram of {\tau}_0')

figure;
histogram(taoNp(1:Npaths,2), 'Normalization', 'probability')
xlabel('{\tau}_1')
ylabel('Probability')
ylim([0 1])
title('Histogram of {\tau}_1')

figure;
histogram(taoNp(1:Npaths,3), 'Normalization', 'probability')
xlabel('Length Stressed Period')
ylabel('Probability')
ylim([0 0.3])
title('Histogram of Length Stressed Period')


%Solution to h
ItNpMax = zeros(Npaths,1);

for i=1:Npaths
    ItNpMax(i,1) = max(ItNp(i,:));
end

ItNpMaxMean = mean(ItNpMax)

figure;
histogram(ItNpMax, 'Normalization', 'probability')
xlabel("Maximum Simulated I's")
ylabel('Probability')
ylim([0 0.3])
title("Histogram of Maximum Simulated I's")


%Solution to i
Prob_OW_D = 0; %Probability that healthcare system will be overwhelmed in Determenistic Case will be zero because Maximum I is Imax

%To calculate Probability that healthcare system will be overwhelmed in Stochastic Case 
Prob_OW_S = 0; 

for i=1:Npaths
    ArrayTemp = find(ItNp(i,:)>1.05*Imaximum);
    if(~isempty(ArrayTemp))
        Prob_OW_S = Prob_OW_S + 1;
    end
end
Prob_OW_S = Prob_OW_S/Npaths


    
