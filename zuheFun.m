% function [MissileLocation,TargetLocation,AM,Los,LosRate,R,dR,MissDistance,count]=zuheFun(TargetType,GuidanceType)
clc;clear;

TargetType=0;
%0  Ŀ��ǻ���
%1  Ŀ���ȼ���
%2  Ŀ��S����
GuidanceType=5;
%0  һ���������
%1  �������Լ����ƫ�ñ���������
%2  ���ڱ�ϵ�����������ĸĽ�����������
%3  �������Ź۲������Ƶ���
%4  �������Լ�������ű������������ڴ�ͳ����������
Flag=1;
% Flag=0; %��ͳ��������

n=2;%nά������
Vt=0;%Ŀ���ٶ�
Vm=180;%�����ٶ�
g=9.8;%�������ٶ�
Vr=zeros(n,1);
V_xiangduiL=zeros(n,1);
V_xiangduiT=zeros(n,1);
Am_dandao=zeros(n,1);%�����ڹ�������ϵ�ļ��ٶȾ���
At_dandao=zeros(n,1);%Ŀ���ڹ�������ϵ�ļ��ٶȾ���
L1=zeros(n,1);
L2=zeros(n,1);
Tp(:,1)=[2000;0];%Ŀ���ʼ����
Mp(:,1)=[0;0];%������ʼ����
sitaMf=-90*pi/180;%�ն�Լ����

Rr=Tp(:,1)-Mp(:,1);%��Ŀ������꣬�ڹ���ϵ
r=sqrt(sum(Rr.^2));%��Ŀ��Ծ���
sitaLOS=atan(Rr(2,1)/sqrt(Rr(1,1)^2));%���߽ǣ���Թ���ϵ��
% sitaM=sitaLOS;
sitaM=0*pi/180;%������ʼ�����
% faiL=atand(-Rr(3,1)/Rr(1,1));

sitaT=140*pi/180;%Ŀ���ʼ�����
% faiT=0;
% ����ϵ
VT(:,1)=[Vt*cos(sitaT);Vt*sin(sitaT)];%Ŀ���ٶ��ڹ�������ϵ��ͶӰ
VM(:,1)=[Vm*cos(sitaM);Vm*sin(sitaM)];%�����ٶ��ڹ�������ϵ��ͶӰ
Vr=VT-VM;%����˶��ٶ��ڹ�������ϵ��ͶӰ
dr=((Rr(1,1)*Vr(1,1))+(Rr(2,1)*Vr(2,1)))/r;%��Ŀ��Ծ���仯��
dsitaLOS=(Vr(2,1)*Rr(1,1)-Vr(1,1)*Rr(2,1))/(r^2);
% %��������һ�ּ��㷽�������ǲ��ԣ������
% dr=Vt*cos(-sitaLOS)-Vm*cos(sitaM-sitaLOS);
% dsitaLOS=(Vt*sin(-sitaLOS)-Vm*sin(sitaM-sitaLOS))/r;

% �����������
T=0.1;
i=1;
wucha = 0.4;
Z=[0;0];
E=[0;0];
u=0;
while(r>10^-3)
    %     if(r<20)
    %         T=0.001;
    %     end
    switch TargetType
        case 0   %Ŀ��ǻ���
            At=0;
        case 1   %Ŀ��ǻ���
            At=30;
        case 2   %Ŀ�����һ���
            if(sin(i*T*pi/5)>0)
                At=-30;
            else
                At=30;
            end
    end
    if(Flag == 1)
        switch GuidanceType
            case 0  %һ���������
                k=4;
                if(r>50)
                    u=k*dsitaLOS*Vm;
                else
                    u=u;
                end
                
            case 1   %�������Լ����ƫ�ñ���������
                k=4;%������ϵ��
                if(r>50)
                    tgo=r/Vm;
                    N=k*Vm;
                    u=N*dsitaLOS+(Vm*(sitaMf-sitaM)-N*(sitaMf-sitaLOS))/tgo;
                else
                    u=u;
                end
            case 2  %���ڱ�ϵ�����������ĸĽ����������ɡ�������
                k=4;%������ϵ��
                if(r>50)
                    u=k*dsitaLOS*Vm/g+k*Vm*cos(sitaM-sitaLOS)*cos(sitaM)/(2*abs(dr));
                else
                    u=u;
                end
            case 3  %�������Ź۲������Ƶ���
                a1=0.5;
                a2=0.25;
                beta01=50;
                beta02=170;
                deta=0.01;
                K=4;
                e=Z(1)-dsitaLOS;
                f0=-2*dr*dsitaLOS/r;
                b0=-1/r;
                dZ=[Z(2)-beta01*e+f0+b0*u;
                    -beta02*fal(e,a1,deta)];
                Z=Z+dZ*T;
                u=K*abs(dr)*dsitaLOS+r*Z(2)/2;
            case 5  %�������Ƕ�Լ��������Ӧ�ն˻�ģ������
                beta=0.2;
                p=5;
                q=3;
                sigma=100;
                alpha = 1;
                gamma = 3;
                k = 100;
                c=5;
                deta= 0.001;
                e1=sitaLOS-sitaMf;
                e2=dsitaLOS;
%                 s=e1+1/beta*e2^(p/q);
                efanshu=abs(e1)+abs(e2);
                s=e1+1/beta*e2^(p/q)+alpha*abs(e1)^gamma*sign(e1);               
                if(r>50)
%                     u=(-2*dr*e2+r*beta*q/p*e2^(2-p/q)+sigma*sign(s))/cos(sitaLOS-sitaM);
                    u=(-2*dr*e2+r*beta*q/p*e2^(2-p/q)*(1+alpha*gamma*abs(e1)^(gamma-1))+sigma*s/((1+c*efanshu)*(abs(s)+deta))+(k+c*efanshu)*s)/(cos(sitaLOS-sitaM));
                else
                    u=u;
                end
                
            case 333 %�������Ź۲������Ƶ���
                alpha1 = 0.5;
                alpha2 = 2/3;
                beta01 = 100;
                beta02 = 300;
                beta03 = 1000;
                k1 = 3;
                k2 = 5;
                tao = 0.45;
                w0 = 60;
                b = -1/(r*tao);
                f0=-3*dr*Z(2)/r+u/(r*tao);
                
                dE=[E(2)-beta01*E(1);
                    E(3)-beta02*sqrt(abs(E(1)))*sign(E(1))-3*dr*E(2)/r;
                    -beta03*((abs(E(1)))^0.25)*sign(E(1))-w0];
                E=E+dE*T;
                dZ=[Z(2)-beta01*E(1);
                    Z(3)-beta02*sqrt(abs(E(1)))*sign(E(1))+b*u+f0;
                    -beta03*((abs(E(1)))^0.25)*sign(E(1))];
                Z=Z+dZ*T;
                
                u=(-k1*sig(Z(3),alpha1)-k2*sig(Z(2),alpha2)-f0-Z(3))/b;
        end
        AT(i)=At;
        AM(i)=u;
        Q(i)=sitaLOS*180/pi;
        dQ(i)=dsitaLOS*180/pi;
        R(i)=r;
        dR(i)=dr;
        Msita(i)=sitaM;
        
        At_dandao(1,i)=AT(i)*cos(sitaT+pi/2);
        At_dandao(2,i)=AT(i)*sin(sitaT+pi/2);
        Am_dandao(:,i)=AM(i)*[cos(sitaM+pi/2);sin(sitaM+pi/2)];
        
        VT(:,i+1)=VT(:,i)+At_dandao(:,i)*T;
        VM(:,i+1)=VM(:,i)+Am_dandao(:,i)*T;
        Tp(:,i+1)=Tp(:,i)+VT(:,i)*T;
        Mp(:,i+1)=Mp(:,i)+VM(:,i)*T;
        Vr=VT(:,i+1)-VM(:,i+1);
        i=i+1;
        %     sitaT=atan(VT(2,i)/VT(1,i));
        sitaM=atan(VM(2,i)/VM(1,i));
        
        
        Rr=Tp(:,i)-Mp(:,i);
        sitaLOS=atan(Rr(2,1)/sqrt(Rr(1,1)^2));%���߽ǣ���Թ���ϵ��
        r=sqrt(sum(Rr.^2));%��Ŀ��Ծ���
        dr=((Rr(1,1)*Vr(1,1))+(Rr(2,1)*Vr(2,1)))/r;%��Ŀ��Ծ���仯��
        dsitaLOS=(Vr(2,1)*Rr(1,1)-Vr(1,1)*Rr(2,1))/(r^2);
        %         dsitaLOS=dsitaLOS+wucha*rand;%������
        
        %��������һ�ּ��㷽���������Ҿ��ò��ԣ������
        %     sitaLOS=atan(Rr(2,1)/sqrt(Rr(1,1)^2));%���߽ǣ���Թ���ϵ��
        %     dr=Vt*cos(-sitaLOS)-Vm*cos(sitaM-sitaLOS);
        %     dsitaLOS=(Vt*sin(-sitaLOS)-Vm*sin(sitaM-sitaLOS))/r;
    else %��ͳ��������
        k=4;
        switch GuidanceType
            case 4  %�������Լ�������ű�������
                if(r>50)
                    %                     dsitaM=k*dsitaLOS+2*abs(dr)*(sitaLOS-sitaMf)/r;
                    dsitaM = k*dsitaLOS;
                else
                    dsitaM=dsitaM;
                end
                %             case 2  %���ڱ�ϵ�����������ĸĽ�����������
                %                 k=4;%������ϵ��
                %                 if(r>50)
                %                     dsitaM=k*dsitaLOS*Vm/g+k*Vm*cos(sitaM-sitaLOS)*cos(sitaM)/(2*abs(dr));
                %                 else
                %                     dsitaM=dsitaM;
                %                 end
                
        end
        
        AT(i)=At;
        AM(i)=Vm*dsitaM;
        Q(i)=sitaLOS*180/pi;
        dQ(i)=dsitaLOS*180/pi;
        R(i)=r;
        dR(i)=dr;
        Msita(i)=sitaM;
        
        At_dandao(1,i)=AT(i)*cos(sitaT+pi/2);
        At_dandao(2,i)=AT(i)*sin(sitaT+pi/2);
        VT(:,i+1)=VT(:,i)+At_dandao(:,i)*T;
        %         VM(:,i+1)=VM(:,i)+Am_dandao(:,i)*T;
        VM(:,i+1)=[Vm*cos(sitaM);Vm*sin(sitaM)];%�����ٶ��ڹ�������ϵ��ͶӰ
        
        Tp(:,i+1)=Tp(:,i)+VT(:,i)*T;
        Mp(:,i+1)=Mp(:,i)+VM(:,i)*T;
        Vr=VT(:,i+1)-VM(:,i+1);
        i=i+1;
        %     sitaT=atan(VT(2,i)/VT(1,i));
        sitaM=sitaM+dsitaM*T;
        
        
        Rr=Tp(:,i)-Mp(:,i);
        sitaLOS=atan(Rr(2,1)/sqrt(Rr(1,1)^2));%���߽ǣ���Թ���ϵ��
        r=sqrt(sum(Rr.^2));%��Ŀ��Ծ���
        dr=((Rr(1,1)*Vr(1,1))+(Rr(2,1)*Vr(2,1)))/r;%��Ŀ��Ծ���仯��
        dsitaLOS=(Vr(2,1)*Rr(1,1)-Vr(1,1)*Rr(2,1))/(r^2);
    end
    if(r<1|i>10000|dr>0)
        r=abs(r-dr*T)
        dr
        break;
    end
    
end
MissileLocation=Mp;
TargetLocation=Tp;
count=i;
MissDistance=r;
Los=Q;
LosRate=dQ;
subplot(2,3,1);
plot(Mp(1,:),Mp(2,:),'green');
hold on
plot(Tp(1,:),Tp(2,:),'red');
title('��Ŀ�˶��켣');
legend('�����˶��켣','Ŀ���˶��켣',2);
subplot(2,3,2);
plot(AM(1:i-1));
title('�������ٶ�');
subplot(2,3,3);
plot(AT(1:i-1));
title('Ŀ����ٶ�');
subplot(2,3,4)
plot(Los(1:i-1));
title('�������');
subplot(2,3,5)
plot(LosRate(1:i-1));
title('������Ǳ仯��');
subplot(2,3,6)
plot(R(1:i-1));
title('��Ŀ����');





