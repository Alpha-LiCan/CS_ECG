function [ ww ] = dwtmtx( N,wtype,wlev )
%DWTMTX Discrete wavelet transform matrix ��ɢС���任����
%   �˺������ɱ任����ww
%   ���� N,wtype,wlev .
%��ϸ��������:
%   N �� ww ��ά��
%   wtype ��С����������
%   wlev �Ƿֽ����
%NOTE: �ź���չģʽ����Ϊ Periodization('per')
[h,g]= wfilters(wtype,'d');         % ���� ��ͨ&��ͨ �˲��� 
L=length(h);                        %�˲�������
h_1 = fliplr(h);                    %���˲���������
g_1 = fliplr(g);
loop_max = log2(N);                 %������
loop_min = double(int8(log2(L)))+1; %��С����
if wlev>loop_max-loop_min+1
    fprintf('\n���棺�ֽ��������\n');
    fprintf('���ֽ����Ϊ %d\n',loop_max-loop_min+1);
    wlev = loop_max-loop_min+1;
end
ww=1;   %Ԥ�������
%%������
for loop = loop_max-wlev+1:loop_max  
    Nii = 2^loop;  
    p1_0 = [h_1 zeros(1,Nii-L)];  %��������
    p2_0 = [g_1 zeros(1,Nii-L)];  
    p1 = zeros(Nii/2,Nii);  
    p2 = zeros(Nii/2,Nii);  
    for ii=1:Nii/2  
        p1(ii,:)=circshift(p1_0',2*(ii-1)+1-(L-1)+L/2-1)';  %ѭ����λ
        p2(ii,:)=circshift(p2_0',2*(ii-1)+1-(L-1)+L/2-1)';  
    end  
    w1=[p1;p2];  %������������
    mm=2^loop_max-length(w1);  
    w=[w1,zeros(length(w1),mm);zeros(mm,length(w1)),eye(mm,mm)];  %���ɵ�λ����
    ww=ww*w;  
    clear p1;clear p2;  
end 
