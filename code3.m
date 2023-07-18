clc        %С������ѡ��
 
clear;    
           
load('ecg1.mat');
x=m;
 
N=1024;
M=512;
k=zeros(52,9);
%52��С����
wav(1)={'haar'}; %���ֽ����Ϊ 9         
wav(2)={'db2'};  %���ֽ����Ϊ 8
wav(3)={'db3'};  %���ֽ����Ϊ 7
wav(4)={'db4'};  %���ֽ����Ϊ 7 
wav(5)={'db5'};  %���ֽ����Ϊ 7
wav(6)={'db6'};  %���ֽ����Ϊ 6
wav(7)={'db7'};  %���ֽ����Ϊ 6
wav(8)={'db8'};  %���ֽ����Ϊ 6
wav(9)={'db9'};  %���ֽ����Ϊ 6
wav(10)={'db10'};%���ֽ����Ϊ 6
wav(11)={'sym2'};%���ֽ����Ϊ 8
wav(12)={'sym3'};%���ֽ����Ϊ 7
wav(13)={'sym4'};%���ֽ����Ϊ 7
wav(14)={'sym5'};%���ֽ����Ϊ 7
wav(15)={'sym6'};%���ֽ����Ϊ 6
wav(16)={'sym7'};%���ֽ����Ϊ 6
wav(17)={'sym8'};%���ֽ����Ϊ 6
wav(18)={'coif1'};%���ֽ����Ϊ 7
wav(19)={'coif2'};%���ֽ����Ϊ 6
wav(20)={'coif3'};%���ֽ����Ϊ 6
wav(21)={'coif4'};%���ֽ����Ϊ 5
wav(22)={'coif5'};%���ֽ����Ϊ 5
wav(23)={'bior1.1'};%���ֽ����Ϊ 9
wav(24)={'bior1.3'};%���ֽ����Ϊ 7
wav(25)={'bior1.5'};%���ֽ����Ϊ 7
wav(26)={'bior2.2'};%���ֽ����Ϊ 7
wav(27)={'bior2.4'};%���ֽ����Ϊ 7
wav(28)={'bior2.6'};%���ֽ����Ϊ 6
wav(29)={'bior2.8'};%���ֽ����Ϊ 6
wav(30)={'bior3.1'};%���ֽ����Ϊ 8
wav(31)={'bior3.3'};%���ֽ����Ϊ 7
wav(32)={'bior3.5'};%���ֽ����Ϊ 6
wav(33)={'bior3.7'};%���ֽ����Ϊ 6
wav(34)={'bior3.9'};%���ֽ����Ϊ 6
wav(35)={'bior4.4'};%���ֽ����Ϊ 7
wav(36)={'bior5.5'};%���ֽ����Ϊ 6
wav(37)={'bior6.8'};%���ֽ����Ϊ 6
wav(38)={'rbio1.1'};%���ֽ����Ϊ 9
wav(39)={'rbio1.3'};%���ֽ����Ϊ 7
wav(40)={'rbio1.5'};%���ֽ����Ϊ 7
wav(41)={'rbio2.2'};%���ֽ����Ϊ 7
wav(42)={'rbio2.4'};%���ֽ����Ϊ 7
wav(43)={'rbio2.6'};%���ֽ����Ϊ 6
wav(44)={'rbio2.8'};%���ֽ����Ϊ 6
wav(45)={'rbio3.1'};%���ֽ����Ϊ 8
wav(46)={'rbio3.3'};%���ֽ����Ϊ 7
wav(47)={'rbio3.5'};%���ֽ����Ϊ 6
wav(48)={'rbio3.7'};%���ֽ����Ϊ 6
wav(49)={'rbio3.9'};%���ֽ����Ϊ 6
wav(50)={'rbio4.4'};%���ֽ����Ϊ 7
wav(51)={'rbio5.5'};%���ֽ����Ϊ 6
wav(52)={'rbio6.8'};%���ֽ����Ϊ 6
tic
for i=1:52       % С����  52�ֱ���
    wtype=wav{i};
    for j=2:10   %�ֽ���� 9�ֱ���
        [ww]=dwtmtx( N,wtype,j);
        Psi=[ww];
        y1 = (Psi*x)';
        threshold=0.05;
        K=length(find(abs(y1)>threshold));     %�ź�x��ϡ���
        k(i,j)=K;
    end
    sprintf('%s%c','ϡ���',num2str(wtype),  '���')
end
save K k;
toc

%%��ͼ
clc;clear all;close all
load('K.mat')
S = ['-ks';'-go';'-m+';'-kd';'-gv';'-b*';'-rx';'-cd';'-r+';'-bv';'-g*';'-md';'-bx';'-cp';'-kh';'-mo'];
figure;
wlev =2:9;
for ii = 1:10
    plot(wlev,k(ii,wlev),S(ii,:));%���x�Ļָ��ź�
    hold on;
end
hold off;
axis([2 10 95 300]);
legend('haar','db2','db3','db4','db5','db6','db7','db8','db9','db10');
xlabel('�ֽ����wlev');
ylabel('ϡ���K');
title('HaarС����DaubechiesС��ϡ���  ͼ��a��');
%
figure;
for ii = 1:7
    plot(wlev,k(ii+10,wlev),S(ii,:));%���x�Ļָ��ź�
    hold on;
end
hold off;
axis([2 10 95 300]);
legend('sym2','sym3','sym4','sym5','sym6','sym7','sym8');
xlabel('�ֽ����wlev');
ylabel('ϡ���K');
title('SymletС��ϡ���');
%
figure;
for ii = 1:5
    plot(wlev,k(ii+17,wlev),S(ii,:));%���x�Ļָ��ź�
    hold on;
end
hold off;
axis([2 10 95 300]);
legend('coif1','coif2','coif3','coif4','coif5');
xlabel('�ֽ����wlev');
ylabel('ϡ���K');
title('CoifletС��ϡ���');
%
figure;
for ii = 1:15
    plot(wlev,k(ii+22,wlev),S(ii,:));%���x�Ļָ��ź�
    hold on;
end
hold off;
axis([2 10 89 300]);
legend('bior1.1','bior1.3','bior1.5','bior2.2','bior2.4','bior2.6','bior2.8','bior3.1','bior3.3','bior3.5','bior3.7','bior3.9','bior4.4','bior5.5','bior6.8');
xlabel('�ֽ����wlev');
ylabel('ϡ���K');
title('BiorthogonalС��ϡ���');
%
figure;
for ii = 1:15
    plot(wlev,k(ii+37,wlev),S(ii,:));%���x�Ļָ��ź�
    hold on;
end
hold off;
axis([2 10 95 300]);
legend('rbio1.1','rbio1.3','rbio1.5','rbio2.2','rbio2.4','rbio2.6','rbio2.8','rbio3.1','rbio3.3','rbio3.5','rbio3.7','rbio3.9','rbio4.4','rbio5.5','rbio6.8');
xlabel('�ֽ����wlev');
ylabel('ϡ���K');
title('ReverseBiorС��ϡ���');
