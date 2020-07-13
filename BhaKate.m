%% 正碰
%kate和bha

clear;
clc;
close all;

%% 物性参数
FuelProperty;

%模型参数
g=9.8;
c1=0.6137;
c2=0.4755;
f0=1.402;
epis=1;

%初始化
R0=1e-3;
theta=90/180*pi;

%定义q得到uin或直接定义uin
QNum=20;
Qmax=40/3600/1000;
Qmin=15/3600/1000;
QSet=linspace(Qmin,Qmax,QNum)';
uin=QSet/R0^2/pi;
results=zeros(QNum,2);

%% 斜碰时定义，正碰时保持phiNum=1
phiNum=1;
phi=linspace(0,pi/1,phiNum)';
r0=R0*sin(theta)./(1-cos(theta)*cos(phi));%输出
Qim=r0.^2*uin*sin(theta)*pi;%输出
u0x=uin*sin(theta)*cos(phi)+uin*cos(theta);
u0y=uin*sin(theta)*sin(phi);
u=(u0x.^2+u0y.^2).^0.5;%输出

for step=1:QNum
    
    %初始化变量
    num=1000000;%最大迭代步数
    us=zeros(num,2);
    r=zeros(num,1);
    us_r=zeros(num,2);
    mark=-1*ones(2,1);
    
    %赋初值
    r(1)=r0;%初半径 %待输入值
    q=Qim(step);%流量 %待输入值
    u0=u(step);%径向平均初速度 %待输入值
    
    us(1,1)=u0/c1;%径向特征初速度
    us(1,2)=u0;
    deltar=1e-6;
    for i=2:num
        r(i)=r(1)+(i-1)*deltar;
    end
    
    %微分方程系数
    a=sig;
    b=2*pi*f0*c1*miu/q;
    c=c2*rho*q/2/pi/c1;
    d=rho*g*q^2/8/pi^2/c1^2;
    
    %迭代求解us
    %bhar
    for i=2:num
        ri=r(i-1);
        usi=us(i-1,1);
        us_r(i,1)=(a*usi-b*ri^2*usi^3+d/ri^2/usi)/(c*usi-a*ri-d/ri/usi^2);
        %                 fprintf('%3.4e %3.4e %3.4e\n',c*usi,a*ri,d/ri/usi^2);
        us(i,1)=us(i-1,1)+deltar*us_r(i,1);
        
        %判断是否发生水跃
        if (us_r(i,1)-us_r(i-1,1))<-2000
            mark(1)=i;
            break
        end
    end
    
    e=g*q/2/pi;
    f=4*pi^2*miu/rho*3/q^2;
    h=1.2;
    k=g*q/2/pi;
    %kate
    for i=2:num
        ri=r(i-1);
        usi=us(i-1,2);
        us_r(i,2)=(e/usi/ri^2-f*ri^2*usi^3)/(h*usi-k/ri/usi^2);
        %         fprintf('%3.4e %3.4e \n',h*usi,k/ri/usi^2);
        us(i,2)=us(i-1,2)+deltar*us_r(i,2);
        
        %判断是否发生水跃
        if (us_r(i,2)-us_r(i-1,2))<-2000
            mark(2)=i;
            break
        end
    end
    
    
    
    %储存结果
    if min(mark)>0
        results(step,:)=[r(mark(1)) r(mark(2))];
    else
        fprintf('第%d步无水跃点\n',step);
    end
end
%% 后处理
Post;
disp("calc end");