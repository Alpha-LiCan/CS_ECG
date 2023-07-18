function [ ww ] = dwtmtx( N,wtype,wlev )
%DWTMTX Discrete wavelet transform matrix 离散小波变换矩阵
%   此函数生成变换矩阵ww
%   参数 N,wtype,wlev .
%详细解释如下:
%   N 是 ww 的维度
%   wtype 是小波函数类型
%   wlev 是分解层数
%NOTE: 信号扩展模式必须为 Periodization('per')
[h,g]= wfilters(wtype,'d');         % 生成 低通&高通 滤波器 
L=length(h);                        %滤波器长度
h_1 = fliplr(h);                    %求滤波器的逆序
g_1 = fliplr(g);
loop_max = log2(N);                 %最大层数
loop_min = double(int8(log2(L)))+1; %最小层数
if wlev>loop_max-loop_min+1
    fprintf('\n警告：分解层数过大\n');
    fprintf('最大分解层数为 %d\n',loop_max-loop_min+1);
    wlev = loop_max-loop_min+1;
end
ww=1;   %预处理矩阵
%%矩阵构造
for loop = loop_max-wlev+1:loop_max  
    Nii = 2^loop;  
    p1_0 = [h_1 zeros(1,Nii-L)];  %构造向量
    p2_0 = [g_1 zeros(1,Nii-L)];  
    p1 = zeros(Nii/2,Nii);  
    p2 = zeros(Nii/2,Nii);  
    for ii=1:Nii/2  
        p1(ii,:)=circshift(p1_0',2*(ii-1)+1-(L-1)+L/2-1)';  %循环移位
        p2(ii,:)=circshift(p2_0',2*(ii-1)+1-(L-1)+L/2-1)';  
    end  
    w1=[p1;p2];  %构造正交矩阵
    mm=2^loop_max-length(w1);  
    w=[w1,zeros(length(w1),mm);zeros(mm,length(w1)),eye(mm,mm)];  %生成单位矩阵
    ww=ww*w;  
    clear p1;clear p2;  
end 
