

clear;close all;clc;
 
%------ SPECIFY DATA ------------------------------------------------------
%------ 指定数据文件 -------------------------------------------------------
PATH=  'D:\university\mit-bih-arrhythmia-database-1.0.0'; % 指定数据的储存路径
HEADERFILE= '100.hea';      % .hea 格式，头文件，可用记事本打开
ATRFILE= '100.atr';         % .atr 格式，属性文件，数据格式为二进制数
DATAFILE='100.dat';         % .dat 格式，ECG 数据
SAMPLES2READ=1024;          % 指定需要读入的样本数
                            % 若.dat文件中存储有两个通道的信号:
                            % 则读入 2*SAMPLES2READ 个数据 
 
%------ LOAD HEADER DATA --------------------------------------------------
%------ 读入头文件数据 -----------------------------------------------------
%
% 示例：用记事本打开的117.hea 文件的数据
%
%      117 2 360 650000
%      117.dat 212 200 11 1024 839 31170 0 MLII
%      117.dat 212 200 11 1024 930 28083 0 V2
%      # 69 M 950 654 x2
%      # None
%
%-------------------------------------------------------------------------
fprintf(1,'\\n$> WORKING ON %s ...\n', HEADERFILE); % 在Matlab命令行窗口提示当前工作状态
% 
% 【注】函数 fprintf 的功能将格式化的数据写入到指定文件中。
% 表达式：count = fprintf(fid,format,A,...)
% 在字符串'format'的控制下，将矩阵A的实数数据进行格式化，并写入到文件对象fid中。该函数返回所写入数据的字节数 count。
% fid 是通过函数 fopen 获得的整型文件标识符。fid=1，表示标准输出（即输出到屏幕显示）；fid=2，表示标准偏差。
%
signalh= fullfile(PATH, HEADERFILE);    % 通过函数 fullfile 获得头文件的完整路径
fid1=fopen(signalh,'r');    % 打开头文件，其标识符为 fid1 ，属性为'r'--“只读”
z= fgetl(fid1);             % 读取头文件的第一行数据，字符串格式
A= sscanf(z, '%*s %d %d %d',[1,3]); % 按照格式 '%*s %d %d %d' 转换数据并存入矩阵 A 中
nosig= A(1);    % 信号通道数目
sfreq=A(2);     % 数据采样频率
clear A;        % 清空矩阵 A ，准备获取下一行数据
for k=1:nosig           % 读取每个通道信号的数据信息
    z= fgetl(fid1);
    A= sscanf(z, '%*s %d %d %d %d %d',[1,5]);
    dformat(k)= A(1);           % 信号格式; 这里只允许为 212 格式
    gain(k)= A(2);              % 每 mV 包含的整数个数
    bitres(k)= A(3);            % 采样精度（位分辨率）
    zerovalue(k)= A(4);         % ECG 信号零点相应的整数值
    firstvalue(k)= A(5);        % 信号的第一个整数值 (用于偏差测试)
end
fclose(fid1);
clear A;
 
%------ LOAD BINARY DATA --------------------------------------------------
%------ 读取 ECG 信号二值数据 ----------------------------------------------
%
if dformat~= [212,212], error('this script does not apply binary formats different to 212.'); end;
signald= fullfile(PATH, DATAFILE);            % 读入 212 格式的 ECG 信号数据
fid2=fopen(signald,'r');
A= fread(fid2, [3, SAMPLES2READ], 'uint8')';  % matrix with 3 rows, each 8 bits long, = 2*12bit
fclose(fid2);
% 通过一系列的移位（bitshift）、位与（bitand）运算，将信号由二值数据转换为十进制数
M2H= bitshift(A(:,2), -4);        %字节向右移四位，即取字节的高四位
M1H= bitand(A(:,2), 15);          %取字节的低四位
PRL=bitshift(bitand(A(:,2),8),9);     % sign-bit   取出字节低四位中最高位，向右移九位
PRR=bitshift(bitand(A(:,2),128),5);   % sign-bit   取出字节高四位中最高位，向右移五位
M( : , 1)= bitshift(M1H,8)+ A(:,1)-PRL;
M( : , 2)= bitshift(M2H,8)+ A(:,3)-PRR;
if M(1,:) ~= firstvalue, error('inconsistency in the first bit values'); end;
switch nosig
case 2
    M( : , 1)= (M( : , 1)- zerovalue(1))/gain(1);
    M( : , 2)= (M( : , 2)- zerovalue(2))/gain(2);
    TIME=(0:(SAMPLES2READ-1))/sfreq;
case 1
    M( : , 1)= (M( : , 1)- zerovalue(1));
    M( : , 2)= (M( : , 2)- zerovalue(1));
    M=M';
    M(1)=[];
    sM=size(M);
    sM=sM(2)+1;
    M(sM)=0;
    M=M';
    M=M/gain(1);
    TIME=(0:2*(SAMPLES2READ)-1)/sfreq;
otherwise  % this case did not appear up to now!
    % here M has to be sorted!!!
    disp('Sorting algorithm for more than 2 signals not programmed yet!');
end
clear A M1H M2H PRR PRL;
fprintf(1,'\\n$> LOADING DATA FINISHED \n');
 
%------ LOAD ATTRIBUTES DATA ----------------------------------------------
atrd= fullfile(PATH, ATRFILE);      % attribute file with annotation data
fid3=fopen(atrd,'r');
A= fread(fid3, [2, inf], 'uint8')';
fclose(fid3);
ATRTIME=[];
ANNOT=[];
sa=size(A);
saa=sa(1);
i=1;
while i<=saa
    annoth=bitshift(A(i,2),-2);
    if annoth==59
        ANNOT=[ANNOT;bitshift(A(i+3,2),-2)];
        ATRTIME=[ATRTIME;A(i+2,1)+bitshift(A(i+2,2),8)+...
                bitshift(A(i+1,1),16)+bitshift(A(i+1,2),24)];
        i=i+3;
    elseif annoth==60
        % nothing to do!
    elseif annoth==61
        % nothing to do!
    elseif annoth==62
        % nothing to do!
    elseif annoth==63
        hilfe=bitshift(bitand(A(i,2),3),8)+A(i,1);
        hilfe=hilfe+mod(hilfe,2);
        i=i+hilfe/2;
    else
        ATRTIME=[ATRTIME;bitshift(bitand(A(i,2),3),8)+A(i,1)];
        ANNOT=[ANNOT;bitshift(A(i,2),-2)];
   end
   i=i+1;
end
ANNOT(length(ANNOT))=[];       % last line = EOF (=0)
ATRTIME(length(ATRTIME))=[];   % last line = EOF
clear A;
ATRTIME= (cumsum(ATRTIME))/sfreq;
ind= find(ATRTIME <= TIME(end));
ATRTIMED= ATRTIME(ind);
ANNOT=round(ANNOT);
ANNOTD= ANNOT(ind);
 
%------ DISPLAY DATA ------------------------------------------------------
figure(1); clf, box on, hold on
plot(TIME, M(:,1),'r');
if nosig==2
    plot(TIME, M(:,2),'b');
end
for k=1:length(ATRTIMED)
    text(ATRTIMED(k),0,num2str(ANNOTD(k)));
end
xlim([TIME(1), TIME(end)]);
xlabel('Time / s'); ylabel('Voltage / mV');
legend('MLII导联','V5导联');
string=['ECG signal ',DATAFILE];
title(string);
figure(1);
fprintf(1,'\\n$> DISPLAYING DATA FINISHED \n');
 
% -------------------------------------------------------------------------
fprintf(1,'\\n$> ALL FINISHED \n');
save ecg100 M;
%--------------------------------------------------------------------------
ecg1=M(:,1);%取一列数据，MLII导联
%------------------------------低通滤波器滤除肌电信号--------------------------
Fs=1500; %采样频率
fp=80;fs=100; %通带截止频率，阻带截止频率
rp=1.4;rs=1.6; %通带、阻带衰减
wp=2*pi*fp;ws=2*pi*fs;
[n,wn]=buttord(wp,ws,rp,rs,'s'); %'s’是确定巴特沃斯模拟滤波器阶次和3dB
%截止模拟频率
[z,P,k]=buttap(n); %设计归一化巴特沃斯模拟低通滤波器，z为极点，p为零点和k为增益
[bp,ap]=zp2tf(z,P,k); %转换为Ha,bp为分子系数，ap为分母系数
[bs,as]=lp2lp(bp,ap,wp); %Ha转换为低通Ha(s)并去归一化，bs为分子系数，as为分母系数
[hs,ws]=freqs(bs,as); %模拟滤波器的幅频响应
[bz,az]=bilinear(bs,as,Fs); %对模拟滤波器双线性变换
[h1,w1]=freqz(bz,az); %数字滤波器的幅频响应
m=filter(bz,az,ecg1);
figure
plot(TIME,ecg1,'b');
hold on;
plot(TIME,m,'r');
xlabel('时间（s）');ylabel('电压（mv）');title('低通滤波前后的波形对比');legend('Original','Recovery');grid;
hold off;
N=1500;
n=0:N-1;
mf=fft(ecg1,N); %进行频谱变换（傅里叶变换）
mag=abs(mf);
f=(0:length(mf)-1)*Fs/length(mf); %进行频率变换
figure
subplot(2,1,1)
plot(f,mag);axis([0,1500,1,50]);grid; %画出频谱图
xlabel('频率(HZ)');ylabel('幅值');title('心电信号频谱图');
mfa=fft(m,N); %进行频谱变换（傅里叶变换）
maga=abs(mfa);
fa=(0:length(mfa)-1)*Fs/length(mfa); %进行频率变换
subplot(2,1,2)
plot(fa,maga);axis([0,1500,1,50]);grid; %画出频谱图
xlabel('频率(HZ)');ylabel('幅值');title('低通滤波后心电信号频谱图');
save ecg1 m
