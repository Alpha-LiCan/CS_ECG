function y_r=resid_paper(y,Phi)
%计算y在Phi上的投影残差

%获取矩阵Phi的行数和列数,M没有用
[M,N]=size(Phi);

%判断矩阵(Phi'*Phi)是否可逆
if(rank(Phi'*Phi)~=N)
    error('矩阵不可逆');
end

y_p=Phi*((Phi'*Phi)\eye(N))*Phi'*y;
y_r=y-y_p;
end
