clear ;clc ;
load('ecg1.mat');%�����ź�
x=m;
N = length(x);%�ź�x�ĳ���
%dct��
ft1=dctmtx(N);
y1=ft1*x;
threshold=0.05;
k=find(abs(y1)>threshold);
K=length(k);%�ź�x��ϡ���
ft1(abs(y1)<threshold)=0;
subplot(311);
plot(y1);
xlabel('��������');
ylabel('��ֵ/mV');
title('DCT��ϡ���K');
%dst��
ft2=dstmtx(N);
y2=ft2*x;
threshold=0.05;
k2=find(abs(y2)>threshold);
K2=length(k2);%�ź�x��ϡ���
ft2(find(abs(y2)<threshold))=0;
subplot(312);
plot(y2);
xlabel('��������');
ylabel('��ֵ/mV');
title('DST��ϡ���K');
%dwt��
wtype = 'db5';
wlev=7;
dwtmode('per');
ww = dwtmtx(N,wtype,wlev);
y3 = (ww*x)';
threshold=0.05;
k3=find(abs(y3)>threshold);
K3=length(k3);%�ź�x��ϡ���
subplot(313);
plot(y3);
xlabel('��������');
ylabel('��ֵ/mV');
title('DWT��ϡ���K');