%% �ж���ʽu_rС��0
%
clear;
clc;
close all;
%����
FuelProperty;
c2=1.2;
g=9.8;
%����
QNum=20;
Qmax=25e-5;
Qmin=1.8e-5;
QSet=linspace(Qmin,Qmax,QNum)';
iter_num=8;%��������
results=zeros(QNum,iter_num);
for step=1:QNum
    
    Q=QSet(step);
    r0=5e-3;
    u0=Q/pi/r0^2;
    h0=r0/2;
    
    %%
    calc_num=1000000;%�����㲽��
    deltar=1e-6;%���㲽��
    r=r0:deltar:(r0+(calc_num-1)*deltar);
    h=zeros(calc_num,iter_num);
    u=zeros(calc_num,iter_num);
    h(1,:)=h0;
    u(1,:)=u0;
    h_r=zeros(calc_num,iter_num);%h��rһ�׵�
    h__r=zeros(calc_num,iter_num);%h��r���׵�
    u_r=zeros(calc_num,iter_num);%u��rһ�׵�
    
    jump_mark=zeros(iter_num,1);
    matrix=zeros(calc_num,iter_num);%��ĸ����
    
    %%
    a=g*Q/2/pi;
    b=4*pi^2*miu/rho*3/Q^2;
    c=2*pi*sig/Q;
    
    for iter=1:iter_num
        if iter==1
            end_mark=calc_num;
        else
            end_mark=jump_mark(iter-1)+1;
        end
        for i=1:end_mark
            denominator=c2*u(i,iter)-a/r(i)/u(i,iter)^2;%��ĸ
            matrix(i,iter)=denominator;
            %�ж�
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
    results(step,:)=jump_mark';%�ռ����
end
%%
results=(results-1).*deltar+r0;
writematrix(results,"MyCalcOut.csv");%���

%% ���
% figure;
% for i=1:iter_num
%     plot(r(1:jump_mark(i))*1000,h(1:jump_mark(i),i)*1000);
%     hold on;
% end
% title('ҺĤ�߶�h�ڲ�ͬr�Ϸֲ�');
% xlabel('r(mm)');
% ylabel('h(mm)');
%
% figure;
% plot(((jump_mark-1)*deltar+r0)*1000,'.-');
% title('ˮԾ�뾶R_j����������仯');
% xlabel('��������');
% ylabel('R_j(mm)');
% ylim([0 inf]);

% figure;
% title('�ٶ�');
% for i=1:iter_num
%     semilogy(r(1:jump_mark(i))*1000,u(1:jump_mark(i),i));
%     hold on;
% end

%Post;
disp("calc end");

