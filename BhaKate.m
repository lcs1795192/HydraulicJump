%% ����
%kate��bha

clear;
clc;
close all;

%% ���Բ���
FuelProperty;

%ģ�Ͳ���
g=9.8;
c1=0.6137;
c2=0.4755;
f0=1.402;
epis=1;

%��ʼ��
R0=1e-3;
theta=90/180*pi;

%����q�õ�uin��ֱ�Ӷ���uin
QNum=20;
Qmax=40/3600/1000;
Qmin=15/3600/1000;
QSet=linspace(Qmin,Qmax,QNum)';
uin=QSet/R0^2/pi;
results=zeros(QNum,2);

%% б��ʱ���壬����ʱ����phiNum=1
phiNum=1;
phi=linspace(0,pi/1,phiNum)';
r0=R0*sin(theta)./(1-cos(theta)*cos(phi));%���
Qim=r0.^2*uin*sin(theta)*pi;%���
u0x=uin*sin(theta)*cos(phi)+uin*cos(theta);
u0y=uin*sin(theta)*sin(phi);
u=(u0x.^2+u0y.^2).^0.5;%���

for step=1:QNum
    
    %��ʼ������
    num=1000000;%����������
    us=zeros(num,2);
    r=zeros(num,1);
    us_r=zeros(num,2);
    mark=-1*ones(2,1);
    
    %����ֵ
    r(1)=r0;%���뾶 %������ֵ
    q=Qim(step);%���� %������ֵ
    u0=u(step);%����ƽ�����ٶ� %������ֵ
    
    us(1,1)=u0/c1;%�����������ٶ�
    us(1,2)=u0;
    deltar=1e-6;
    for i=2:num
        r(i)=r(1)+(i-1)*deltar;
    end
    
    %΢�ַ���ϵ��
    a=sig;
    b=2*pi*f0*c1*miu/q;
    c=c2*rho*q/2/pi/c1;
    d=rho*g*q^2/8/pi^2/c1^2;
    
    %�������us
    %bhar
    for i=2:num
        ri=r(i-1);
        usi=us(i-1,1);
        us_r(i,1)=(a*usi-b*ri^2*usi^3+d/ri^2/usi)/(c*usi-a*ri-d/ri/usi^2);
        %                 fprintf('%3.4e %3.4e %3.4e\n',c*usi,a*ri,d/ri/usi^2);
        us(i,1)=us(i-1,1)+deltar*us_r(i,1);
        
        %�ж��Ƿ���ˮԾ
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
        
        %�ж��Ƿ���ˮԾ
        if (us_r(i,2)-us_r(i-1,2))<-2000
            mark(2)=i;
            break
        end
    end
    
    
    
    %������
    if min(mark)>0
        results(step,:)=[r(mark(1)) r(mark(2))];
    else
        fprintf('��%d����ˮԾ��\n',step);
    end
end
%% ����
Post;
disp("calc end");