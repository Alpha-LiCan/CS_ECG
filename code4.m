clc %%运行得到完整的数据
clear;
load('ecg1.mat');
x=m;
N=1024;
% M=922; % 测量数(M>=K*log(N/K),
CR(1)={'Sheet10'};
CR(2)={'Sheet15'};
CR(3)={'Sheet20'};
CR(4)={'Sheet25'};
CR(5)={'Sheet30'};
CR(6)={'Sheet35'};
CR(7)={'Sheet40'};
CR(8)={'Sheet45'};
CR(9)={'Sheet50'};
CR(10)={'Sheet55'};
CR(11)={'Sheet60'};
CR(12)={'Sheet65'};
CR(13)={'Sheet70'};
wav(1)={'db5'}; %最大分解层数为 7
wav(2)={'sym4'}; %最大分解层数为 7
wav(3)={'coif2'}; %最大分解层数为 6
wav(4)={'bior4.4'}; %最大分解层数为 7
wav(5)={'rbio1.5'}; %最大分解层数为 7
dwtmode('per');
A=zeros(5,10);
for c=1:13 % 压缩率 13种遍历
    switch c
        case 1
            M=922; % 10%
        case 2
            M=870; % 15%
        case 3
            M=819; % 20%
        case 4
            M=768; % 25%
        case 5
            M=719; % 30 %
        case 6
            M=666; % 35%
        case 7
            M=614; % 40%
        case 8
            M=563; % 45 %
        case 9
            M=512; % 50%
        case 10
            M=461; % 55%
        case 11
            M=410; % 60%
        case 12
            M=358; % 65%
        case 13
            M=307; % 70%
    end
    for k= 1:5 % 重构算法 5种遍历
        for j=1:1 % 小波基 5种遍历
            PRDsum=0;
            i=0;
            tic
            while i<10 %i<10; %测试次数 10
                Phi=BernoulliMtx( M,N )/sqrt(M);
                y=Phi*x;
                wtype=wav{j};
                [ww]=dwtmtx( N,wtype,6);
                Psi=ww;
                y1 = (Psi*x)';
                threshold=0.05;
                K=length(find(abs(y1)>threshold)); %信号x的稀疏度
                T=Phi*Psi';
                switch k %重构算法选择 K决定
                    case 1
                        hat_s=CS_OMP(y,T,K);
                        hat_x=real(Psi.*hat_s);
                    case 2
                        hat_s=CS_CoSaMP(y,T,K);
                        hat_x=real(Psi'*hat_s);
                    case 3
                        hat_s=CS_SP(y,T,K);
                        hat_x=real(Psi.*hat_s);
                    case 4
                        hat_s=CS_StOMP(y,T,K);
                        hat_x=real(Psi'*hat_s);
                    case 5
                        hat_s=CS_gOMP(y,T,K);
                        hat_x=real(Psi'*hat_s);
                end
                PRD=norm(x-hat_x)/norm(x)*100;
                PRDsum=PRDsum+PRD;
                i=i+1;
            end % 运行10次 结束
            A(j,2*k-1)=PRDsum/i ;
            A(j,2*k)=toc/i ;
            sprintf('%s%c%s%c%s%c','Sheet' ,num2str(c), '重构算法',num2str(k), '稀疏基',num2str(wtype), '完成')
        end %小波遍历结束
        sprintf('%s%c%s%c','Sheet' ,num2str(c), '重构算法',num2str(k), '完成')
    end % 重构算法遍历结束
    sprintf('%s%c','Sheet', num2str(c), '完成')
    tabel=CR{c};
    xlswrite('test.xlsx',A,tabel,'C4');
    sprintf('Sheet',c, '已存入表 ')
    % 13种压缩率 遍历结束
end
toc
sprintf('All Complete')