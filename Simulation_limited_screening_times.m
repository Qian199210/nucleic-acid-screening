function  Parest_SXNOA  
clear all
clc

% Data_daily=xlsread('Xian','B2:B33');
% data_yuanshi=xlsread('Xian','C2:C33');  %accumulative confirmed
% Data1_yuanshi=xlsread('Xian','D2:D33');%daily confirmed from I
% Data2_yuanshi=xlsread('Xian','E2:E33');%daily confirmed from E_q
% Data3=xlsread('Xian','F2:F53');%confirmed from screening
% size(Data3)
% 
% Num1=1;Num2=32;
% data1=zeros(5,length(data_yuanshi(Num1:Num2)));
% data1(1,:)=0:length(data_yuanshi(Num1:Num2))-1;
% data1(2,:)=data_yuanshi(Num1:Num2);
% data1(3,:)=Data1_yuanshi(Num1:Num2);
% data1(4,:)=Data2_yuanshi(Num1:Num2);
% data1(5,:)=Data3(Num1:Num2);

sigma=1/3;lambda=1/14;  
Sq0=0; Eq0=0; S0=1.29*10^7; R0=0; H0=1; HA0=1; T_switch=13; T_step=1; Run_step=0.01;

par2=[sigma lambda S0 Sq0 Eq0 R0 H0 HA0 T_switch T_step Run_step];

load 'Par_CI'

par1=mean(Par_CI);


beta=par1(1);
r1=par1(2); 
c_0=par1(3);

q_0=par1(5);    
r2=par1(6);
par1(7)=0.3;
q_m=par1(7);
gam_I=par1(8)
delta_I0=par1(9);        
delta_q0=par1(10);
gam_H=par1(11)
I0=par1(12);
E0=par1(13);
delta_I1=par1(14);
delta_q1=par1(15);

sigma=par2(1);
lambda=par2(2);
S0=par2(3); 
Sq0=par2(4);
Eq0=par2(5);
R0=par2(6); 
H0=par2(7); 
HA0=par2(8);
T_switch=par2(9);
T_step=par2(10);
Run_step=par2(11);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%情形1%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par1(4)=c_0*0.9;
c_b=par1(4);

T_final=365;
tspan=0:0.1:T_switch;
x0=[S0 E0 I0 Sq0 Eq0 H0 R0 HA0];
[T,x]=ode45(@kuzode,tspan,x0,[],par1, par2);

figure(1)
subplot(221)
plot(T,x(:,3))
hold on
subplot(222)
plot(T,x(:,6));
hold on

x0_ini=x(end,:);
h=1;
%T_int=1;
TT=T_switch;

CC=[1 0 0;
    0 1 0;
    0 0 1;
    1 0 1];

T_ini=[14 30 45];
qs=0.6;
T_seq=1:0.1:60;

for i=1:length(T_seq)
    i
x0=x0_ini;
h=1;
T_int=T_seq(i);
TT=T_switch;
AA=[];
while TT<T_final%&h<=10
    tspan=[T_switch+(h-1)*T_int:0.1:T_switch+h*T_int];
    [T,x]=ode45(@kuzode,tspan,x0,[],par1, par2);
% figure(1)
% subplot(221)
% plot(T,x(:,3),'color',CC(i,:))
% hold on
% subplot(222)
% plot(T,x(:,6),'color',CC(i,:));
% hold on
    x0=x(end,:);
    x0(6)=x0(6)+qs*x0(3);
    x0(8)=x0(8)+qs*x0(3);  
    x0(3)=(1-qs)*x0(3);
    h=h+1;
    TT=T_switch+h*T_int;
    AA=[AA; x(:,6)];  
end
   Max(i)=max(AA);
% tspan=[T_switch+(h-1)*T_int T_final];
% [T,x]=ode45(@kuzode,tspan,x0,[],par1, par2);
% figure(1)
% subplot(221)
% plot(T,x(:,3),'color',CC(i,:))
% hold on
% subplot(222)
% plot(T,x(:,6),'color',CC(i,:));
% hold on
end

figure(2)
plot(T_seq,Max)

% %%%%%%%%%%%%%%%%%%%%%%%%情形2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% par1(4)=c_0*0.6;
% c_b=par1(4);
% 
% T_final=365;
% tspan=0:0.1:T_switch;
% x0=[S0 E0 I0 Sq0 Eq0 H0 R0 HA0];
% [T,x]=ode45(@kuzode,tspan,x0,[],par1, par2);
% 
% figure(1)
% subplot(223)
% plot(T,x(:,3))
% hold on
% xlabel('Days')
% ylabel('I(t)')
% title('(C)')
% subplot(224)
% plot(T,x(:,6));
% hold on
% xlabel('Days')
% ylabel('H(t)')
% title('(D)')
% x0_ini=x(end,:);
% h=1;
% %T_int=1;
% TT=T_switch;
% 
% CC=[1 0 0;
%     0 1 0;
%     0 0 1;
%     1 0 1];
% 
% qs_ini=[0 0.2 0.7 0.9];
% 
% for i=1:3
%     i
% qs=qs_ini(i)
% x0=x0_ini;
% h=1;
% T_int=7;
% TT=T_switch;
% while TT<T_final%&h<=7
%     tspan=[T_switch+(h-1)*T_int:0.1:T_switch+h*T_int];
%     [T,x]=ode45(@kuzode,tspan,x0,[],par1, par2);
% figure(1)
% subplot(223)
% plot(T,x(:,3),'color',CC(i,:))
% hold on
% subplot(224)
% plot(T,x(:,6),'color',CC(i,:));
% hold on
%     x0=x(end,:);
%     x0(6)=x0(6)+qs*x0(3);
%     x0(8)=x0(8)+qs*x0(3);  
%     x0(3)=(1-qs)*x0(3);
%     h=h+1;
%     TT=T_switch+h*T_int;
% end
% 
% end




% tspan=TT:0.1:365;
% [T,x]=ode45(@kuzode,tspan,x0,[],par1, par2);
% plot(T,x(:,6))










% [T,X,X_fit]=systemsoln(par1,y0,par2,T_final,I_pulse);
% 
% 
% delta_q=delta_q0.*((0:Run_step:T_final)<=T_switch)+delta_q1.*((0:Run_step:T_final)>T_switch);
% delta_I=delta_I0.*((0:Run_step:T_final)<=T_switch)+delta_I1.*((0:Run_step:T_final)>T_switch);
% 
% 
% Fit1_Xian=delta_I'.*X(:,3);
% Fit2_Xian=delta_q'.*X(:,5);
% Fit3_Xian=X(:,8);
% % Sol2(:,i)=X(:,2);
% % Sol3(:,i)=X(:,3);
% % Sol5(:,i)=X(:,5);
% % Sol6(:,i)=X(:,6);
% 
% Daily_confirm=[Fit3_Xian(1,1); diff(Fit3_Xian(:,1))];
% 
% 
% xx=0:Run_step:T_final;
% figure(1)
% subplot(322)
% 
% plot(xx,Fit1_Xian,'color',[240 59 32]/255,'LineWidth',1)
% %plot(T,delta_I'.*X(:,3),'r-')
% hold on
% ylabel('Opportunistic confirmed')
% title('(B)')
% % set(gca,'XTick',[0:6:48])
% % set(gca,'xticklabel',{'12/9','12/15','12/21','12/27','2022/1/2','1/8','1/14','1/20','1/26'});
% %axis([0 32 0 40])
% 
% 
% subplot(324)
% plot(xx,Fit2_Xian,'color',[240 59 32]/255,'LineWidth',1)
% hold on
% % %plot(T,delta_q'.*X(:,5),'r-')
% % hold on
% % plot(data1(1,:),data1(4,:),'ko','Markersize',4)
% % hold on
% % ylabel('Confirmed from E_q')
% % title('(D)')
% % set(gca,'XTick',[0:6:48])
% % set(gca,'xticklabel',{'12/9','12/15','12/21','12/27','2022/1/2','1/8','1/14','1/20','1/26'});
% %axis([0 32 0 150])
% 
% subplot(326)
% 
% plot(xx,Fit3_Xian,'color',[240 59 32]/255,'LineWidth',1)%plot(T,X(:,8),'r-')
% hold on
% 
% % ylabel('Accumulative confirmed cases')
% % xlabel('Days')
% % title('(F)')
% % set(gca,'XTick',[0:6:48])
% % set(gca,'xticklabel',{'12/9','12/15','12/21','12/27','2022/1/2','1/8','1/14','1/20','1/26'});
% %axis([0 32 0 2500])
% 
% %sgtitle('Xian')
% figure(2)
% %ave=mean(Daily_confirm,2);
% % zhixingqujian = [quantile(Daily_confirm',0.025)' , quantile(Daily_confirm',0.975)'];
% % fill([xx,fliplr(xx)],[zhixingqujian(:,1)' ,fliplr(zhixingqujian(:,2)')],[255 237 160]/255,'edgecolor',[255 237 160]/255,'facealpha',0.7,'edgealpha',0.7)
% % hold on 
% plot(xx,X(:,3),'color',[240 59 32]/255,'LineWidth',1)


function xdot=kuzode(tn,x,par1,par2)
beta=par1(1);
r1=par1(2); 
c_0=par1(3); 
c_b=par1(4);
q_0=par1(5);     
r2=par1(6);
q_m=par1(7);
gam_I=par1(8);
delta_I0=par1(9);       
delta_q0=par1(10);
gam_H=par1(11);
I0=par1(12);
E0=par1(13);
delta_I1=par1(14);
delta_q1=par1(15);

sigma=par2(1);
lambda=par2(2);
S0=par2(3); 
Sq0=par2(4);
Eq0=par2(5);
R0=par2(6); 
H0=par2(7); 
HA0=par2(8);
T_switch=par2(9);
Run_step=par2(11);


if tn<=T_switch
    q=q_0;
    c=c_0;
delta_I=delta_I0;
delta_q=delta_q0;
else
%c=(c_0-c_b).*exp(-r1.*(tn-T_switch))+c_b;
c=c_b;
%q=q_m;
delta_I=delta_I1;
delta_q=delta_q1;
q=(q_0-q_m).*exp(-r2.*(tn-T_switch))+q_m;
end


x1=x(1);x2=x(2);x3=x(3);x4=x(4);x5=x(5);x6=x(6);x7=x(7);

N=x1+x2+x3+x4+x5+x6+x7;

xdot=[-beta*c*x1*x3/N-(1-beta)*c*q*x1*(x3)/N+lambda*x4;
      beta*c*(1-q)*x1*x3/N-sigma*x2;
      sigma*x2-delta_I*x3-gam_I*x3;    
      (1-beta)*c*q*x1*x3/N-lambda*x4; 
      beta*c*q*x1*x3/N-delta_q*x5;
      delta_I*x3+delta_q*x5-gam_H*x6;
      gam_I*x3+gam_H*x6;
      delta_q*x5+delta_I*x3];%��
  
  



     
function [T_plot,X,X_fit]=systemsoln(par1,x0,par2,T_final,I_pulse)

T_step=par2(10);
Run_step=par2(11);


tspan=[0:Run_step:T_step];
h=1;
X=[];
T_plot=[];
TT=0;
X_fit=x0;

while TT<T_final

[T,x]=ode45(@kuzode,tspan,x0,[],par1, par2);
    
   if TT<=13| TT/5==0
       qs=0;
   else
       qs=0.3;
   end

    x0=x(end,:);
    x0(6)=x0(6)+qs*x0(3);
    x0(8)=x0(8)+qs*x0(3);  
    x0(3)=(1-qs)*x0(3);
    tspan=[h*T_step:Run_step:(h+1)*T_step];
    h=h+1;
    TT=T(end);
    X=[X; x(1:end-1,:)];
    T_plot=[T_plot; T(1:end-1)];
    X_fit=[X_fit; x0];
    
end

X=[X; x(end,:)];
T_plot=[T_plot; T(end)];





