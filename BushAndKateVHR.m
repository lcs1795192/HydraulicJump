%% 模板：判定方式u_r小于0
%改u_r,h_r,h__r
clear;
clc;
close all;
%常数
c2=1.2;
FuelProperty;
g=9.8;
%条件
qNum=12;
qSet=linspace(2.5e-5,25e-5,qNum)';
iter_num=8;%迭代次数
results=zeros(qNum,iter_num);
for k=1:qNum
    
%q=5e-5;
q=qSet(k);%转换单位
r0=5e-3;
u0=q/pi/r0^2;
h0=r0/2;

%%
calc_num=8000;%最大计算步数
deltar=0.00005;%计算步长
r=r0:deltar:(r0+(calc_num-1)*deltar);
h=zeros(calc_num,iter_num);
u=zeros(calc_num,iter_num);
h(1,:)=h0;
u(1,:)=u0;
h_r=zeros(calc_num,iter_num);%h对r一阶导
h__r=zeros(calc_num,iter_num);%h对r二阶导
u_r=zeros(calc_num,iter_num);%u对r一阶导

jump_mark=zeros(iter_num,1);
matrix=zeros(calc_num,iter_num);%分母矩阵

%%
a=g*q/2/pi;
b=4*pi^2*miu/rho*3/q^2;
s=2*pi*sig/q;

for j=1:iter_num
    if j==1
        end_mark=calc_num;
    else
        end_mark=jump_mark(j-1)+1;
    end
    for i=1:end_mark
        denominator=c2*u(i,j)-a/r(i)/u(i,j)^2;%分母
        matrix(i,j)=denominator;
        %判定
        if denominator<0
            jump_mark(j)=i-1;
            break
        end
        if j==1
            u_r(i,j)=(a/r(i)^2/u(i,j)-b*r(i)^2*u(i,j)^3)/denominator;
        else
            u_r(i,j)=(a/r(i)^2/u(i,j)-b*r(i)^2*u(i,j)^3-...
                s*u(i,j)*h_r(i,j-1)*(r(i)*h__r(i,j-1)+h_r(i,j-1)+h_r(i,j-1)^3)/(1+h_r(i,j-1)^2)^1.5)...
                /denominator;
        end
        u(i+1,j)=u(i,j)+deltar*u_r(i,j);
        h_r(i,j)=q/2/pi*(-1/r(i)^2/u(i,j)-1/r(i)/u(i,j)^2*u_r(i,1));
        h(i,j)=q/2/pi/r(i)/u(i,j);
    end
    
    for i=2:jump_mark(j)-1
        h__r(i,j)=(h(i+1,j)+h(i-1,j)-2*h(i,j))/deltar^2;
    end
    h__r(1,j)=(h_r(2,j)-h_r(1,j))/deltar;
    %h__r(mark(j))=(h_r(mark(j))-h_r(mark(j)-1))/deltar;
end
results(k,:)=jump_mark';%收集结果
end
%%
results=((results-1).*deltar+r0)*1000;%转换单位为mm
writematrix(results,"output.txt");%输出

%% 输出
% figure;
% for i=1:iter_num
%     plot(r(1:jump_mark(i))*1000,h(1:jump_mark(i),i)*1000);
%     hold on;
% end
% title('液膜高度h在不同r上分布');
% xlabel('r(mm)');
% ylabel('h(mm)');
% 
% figure;
% plot(((jump_mark-1)*deltar+r0)*1000,'.-');
% title('水跃半径R_j随迭代次数变化');
% xlabel('迭代次数');
% ylabel('R_j(mm)');
% ylim([0 inf]);

% figure;
% title('速度');
% for i=1:iter_num
%     semilogy(r(1:jump_mark(i))*1000,u(1:jump_mark(i),i));
%     hold on;
% end
disp("calc end");

