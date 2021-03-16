%% 判定方式u_r小于0
%
clear;
clc;
close all;
%常数
g=9.8;
FuelProperty;
c2=1.2;
f_0=3;
% type=5;
% switch type
%     case 2
%         c2=1.2;
%         f_0=3;
%     case 3
%         c2=1.243429;
%         f_0=2.4;
%     case 4
%         c2=1.268861;
%         f_0=20/9;
%     case 5
%         c2=1.284919;
%         f_0=15/7;
% end
%条件
QNum=5;
Qmax=45/1000/3600;
Qmin=15/1000/3600;
QSet=linspace(Qmin,Qmax,QNum)';
iter_num=8;%迭代次数
results=zeros(QNum,iter_num);
for step=1:QNum
    
    Q=QSet(step);
    r0=5e-3;
    u0=Q/pi/r0^2;
    h0=r0/2;
    
    %%
    calc_num=1000000;%最大计算步数
    deltar=0.5e-5;%计算步长
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
    a=g*Q/2/pi;
    b=4*pi^2*miu/rho*f_0/Q^2;
    c=2*pi*sig/Q;
    
    for iter=1:iter_num
        if iter==1
            end_mark=calc_num;
        else
            end_mark=jump_mark(iter-1)+1;
        end
        for i=1:end_mark
            denominator=c2*u(i,iter)-a/r(i)/u(i,iter)^2;%分母
            matrix(i,iter)=denominator;
            %判定
            if denominator<0
                jump_mark(iter)=i-1;
                break
            end
            if iter==1
                u_r(i,iter)=(a/r(i)^2/u(i,iter)-b*r(i)^2*u(i,iter)^3)/denominator;
            else
                u_r(i,iter)=(a/r(i)^2/u(i,iter)-b*r(i)^2*u(i,iter)^3-...
                    c*u(i,iter)*h_r(i,iter-1)*(r(i)*h__r(i,iter-1)+h_r(i,iter-1)+h_r(i,iter-1)^3)/(1+h_r(i,iter-1)^2)^1.5)...
                    /denominator;
            end
            u(i+1,iter)=u(i,iter)+deltar*u_r(i,iter);
            h_r(i,iter)=Q/2/pi*(-1/r(i)^2/u(i,iter)-1/r(i)/u(i,iter)^2*u_r(i,1));
            h(i,iter)=Q/2/pi/r(i)/u(i,iter);
        end
        
        for i=2:jump_mark(iter)-1
            h__r(i,iter)=(h(i+1,iter)+h(i-1,iter)-2*h(i,iter))/deltar^2;
        end
        h__r(1,iter)=(h_r(2,iter)-h_r(1,iter))/deltar;
        %h__r(mark(j))=(h_r(mark(j))-h_r(mark(j)-1))/deltar;
    end
    results(step,:)=jump_mark';%收集结果
    
    figure(1);
    hold on;
    plot(r(1:jump_mark(3))*1000,h(1:jump_mark(3),3)*1000);
    writematrix([r(1:jump_mark(3))'*1000,h(1:jump_mark(3),3)*1000],"filmHeight"+step+".csv");%输出
    title('液膜厚度h在不同r上分布');
    xlabel('r(mm)');
    ylabel('h(mm)');
end
%%
results=(results-1).*deltar+r0;
writematrix(results,"MyCalcOut.csv");%输出

%% 输出




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

%%
% Post;
disp("calc end");

