clear ;clc ;
load('ecg1.mat');%导入信号
x=m;
N = length(x);%信号x的长度
%dct基
ft1=dctmtx(N);
y1=ft1*x;
threshold=0.05;
k=find(abs(y1)>threshold);
K=length(k);%信号x的稀疏度
ft1(abs(y1)<threshold)=0;
subplot(311);
plot(y1);
xlabel('采样点数');
ylabel('幅值/mV');
title('DCT域稀疏度K');
%dst基
ft2=dstmtx(N);
y2=ft2*x;
threshold=0.05;
k2=find(abs(y2)>threshold);
K2=length(k2);%信号x的稀疏度
ft2(find(abs(y2)<threshold))=0;
subplot(312);
plot(y2);
xlabel('采样点数');
ylabel('幅值/mV');
title('DST域稀疏度K');
%dwt基
wtype = 'db5';
wlev=7;
dwtmode('per');
ww = dwtmtx(N,wtype,wlev);
y3 = (ww*x)';
threshold=0.05;
k3=find(abs(y3)>threshold);
K3=length(k3);%信号x的稀疏度
subplot(313);
plot(y3);
xlabel('采样点数');
ylabel('幅值/mV');
title('DWT域稀疏度K');