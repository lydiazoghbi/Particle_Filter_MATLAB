clc, clear all, close all 
% =============== Setup
WorkName = 'Battery';
TimeUnit= 'Cycles';
dt=1;

% =============== Parameters
n=1001;
signiLevel=5;
thres= 0.8;
finalPoint = 200;
samplingMethod = 4; % Multinomial, Residual, Stratified, Systematic (For Now)
methods = {'Multinomial','Residual','Stratified','Systematic'};
resampFrequency = 1;

% =============== Data Input / Normalization
inputData = csvread('cx2_37.csv');
data = inputData./inputData(1);
data = smooth(data,.1,'rloess');
% A = 1:.1:100;
% led = 1-A/100;
measuData=data(1:finalPoint);
v = .1;
% =============== Model Parameters
ParamName=['x';'A';'B';'s'];
initDisPar=[.95 1.05; .95 1.05; -.005 0; .075 .2];
%        a = 0.004871  (-0.009499, -0.0002427)
%        b =    0.002656  (0.002024, 0.003287)
%        c =       1.313  (1.309, 1.316)
%        d =  -0.0002425  (-0.0002615, -0.0002236)
% =============== Initial Distribution of Parameters
p=size(ParamName,1);
ParamResul=[];
for j=1:p
    param(j,:)=unifrnd(initDisPar(j,1),initDisPar(j,2),...
        1,n);
    ParamResul=[ParamResul; ParamName(j,:) 'Resul'];
    eval([ParamResul(j,:) '=param(j,:);']);
end

% =============== Positive Trend or Negative Trend?
k1=length(measuData);
k=1;
if measuData(end)-measuData(1)<0
    cofec=-1;
else
    cofec =1;
end
a = 1;
% =============== Main Loop
while min(eval([ParamResul(1,:) '(k,:)'])*cofec)<thres*cofec
  k=k+1;
  paramPredi=param;
  for j=1:p
      eval([ParamName(j,:) '=paramPredi(j,:);']);
  end
  paramPredi(1,:)= A.*exp(B.*k)+v*normrnd(0,.333);
  
  if k<=k1 
      if a == resampFrequency;
          likel=normpdf(paramPredi(1,:),measuData(k),...
              paramPredi(end,:));
          loca = resampleMethod(likel,n,3);
          for i=1:n
              param(:,i) = paramPredi(:,loca(i));
          end
          a = 1;
      else
          a = a+1;
      end
  else
      param=paramPredi;
  end
  
  for j=1:p
      eval([ParamResul(j,:) '(k,:)=param(j,:);']);
  end
  
  if k>k1
      eval([ParamResul(1,:) '(k,:)=normrnd(param(1,:),param(end,:));']);
  end
  
end
% ================== Determine RUL
time=[0:dt:dt*k-1]';
% RUL = zeros(n,1);
perceValue=[50 signiLevel 100-signiLevel];
% for i=1:n
%     loca=find(eval([ParamResul(1,:) '(:,i)'])*cofec>=thres*cofec);
%     if k <= k
%         RUL(i)=time(loca(1))-k1;
%     else
%         RUL(i) = 0;
%     end
% end
% % =================== Histogram Plot
% RULPerce=prctile(RUL', perceValue);
% figure
% set(gca,'fontsize',14);
% hist(RUL,30);
% % xlim([min(RUL) max(RUL)]);
% xlabel(['RUL' '(' TimeUnit ')']);
methodTitle = methods(samplingMethod);
titleName=[methodTitle,' Resampling'];
% title(titleName)
% ylabel('Frequency')
% fprintf('\n # Percentiles of RUL at %g hours \n', k1)
% fprintf('\n %gprct: %g, median: %g, %gprct: %g \n', perceValue(2),...
%     RULPerce(2), RULPerce(1),perceValue(3), RULPerce(3))
% % Name=[WorkName 'at' num2str(time(k1)) '.mat'];
%%
% ==================== Plot the Filter Over Time
figure(2)
plot(repmat(time,1,3),prctile(xResul',perceValue)','o','MarkerSize',3) % Percentiles
hold on
t=0:1:1*length(data)-1;
plot(t',data,'k'); % Actual Data
plot([0,length(data)],[thres,thres],'k:') % Threshold Line
plot([finalPoint,finalPoint],[.6,1.2],'c:') % Provided Data Cutoff
legend('Median Value','5 Percentile','95 Percentile',...
    'True Value','Threshold Value','Training Cutoff',...
    'Location','Best')
title(titleName)
xlabel(TimeUnit)
xlim([0 k])
ylabel('Percentage of Initial Value')
ylim([0.4,1.1])
%% Individual Particle Instances
figure(4)
hold on
for i = 1:50
    val = length(unique(xResul(i,:)));
    scatter(i,val,'k')
end
%% RUL Estimation at each particle at each iteration 
TUL = zeros(k1+50,n);
for i=1:k1+50
    point = [AResul(i,:);BResul(i,:)];
    TUL(i,:) = solverFunSingle(point,thres);
end
medians = median(TUL,2);
RUL = medians - [1:k1+50]';
plot(RUL)
% hold  on
% plot(data,'k')