clear ; close all;
load('ecg1.mat');
y=m;       %载入原始信号
tic
n=length(y);
Fs=360;    %采样频率
Dt=1/Fs;   %时域采样间隔
t=0:Dt:(n-1)/Fs;
f=(0:length(y)-1)*Fs/length(y);%频率变换
T=2;
Dt0=T*Dt;
Fs0=1/Dt0;
t0=0:Dt0:(n-1)/Fs0;
 
for i=1:n
    if mod(i,T)==0        %作除求余
        q=1;
        y1(i)=y(i)*q;     %隔一取一
    else
        q=0;
         y1(i)=y(i)*q;
    end
end
n1=1:n;
nTs=n1*Dt0;
y7 = y1 * sinc( Fs0 * ( ones(length(n1),1) * t - nTs' * ones(1,length(t) ) ));%内插法重构信号
%%
Fw=fft(y,n);        %进行频谱变换（傅里叶变换）
y2=abs(Fw)/n;         
Fw1=fft(y1,n);
y4=abs(Fw1)/n;
Fw2=fft(y7,n);
y6=abs(Fw2)/n;
toc
%%
figure;
subplot(2,2,1);
plot(y);
xlabel(['采样点' 10 '(a)']);
ylabel('电压（mV）');
grid;title('时域连续信号图');
subplot(2,2,2);
plot(f,y2);
xlabel(['频率（Hz）' 10 '(b)']);
ylabel('Xa(jW)');
title('时域信号频谱图');
grid;axis([0,400,0,0.10]);
subplot(2,2,3);
plot(y7);
xlabel(['采样点' 10 '(c)']);
ylabel('电压（mV）');
grid;title('时域内插函数重构信号波形图');
subplot(2,2,4);
plot(f,y6);
xlabel(['频率（Hz）' 10 '(d)']);
ylabel('Xa(jW)');
grid;axis([0,400,0,0.10]);
title('时域内插函数重构信号频谱图');
%%
figure;
subplot(2,2,1);
plot(y);
xlabel(['采样点' 10 '(a)']);
ylabel('电压（mV）');
grid;title('时域连续信号');
subplot(2,2,2);
plot(f,y2);
xlabel(['频率（Hz）' 10 '(b)']);
ylabel('Xa(jW)');
grid;axis([0,400,0,0.10]);
title('时域信号频谱图');
subplot(2,2,3);
stem(y1);
xlabel(['采样点' 10 '(c)']);
ylabel('电压（mV）');
grid;title('采样离散信号');
%plot(y3);
subplot(2,2,4);
plot(f,y4);
xlabel(['频率（Hz）' 10 '(d)']);
ylabel('Xa(jW)');
grid;axis([0,400,0,0.10]);
title('离散信号频谱图');
%snr2=snr(y,y2);
%str='snr_传统方法=';
%str=[str,num2str(snr2)];
%disp(str);
PRD=norm(y-y7)/norm(y)*100 

