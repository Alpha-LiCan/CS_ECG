function y_r=resid_paper(y,Phi)
%����y��Phi�ϵ�ͶӰ�в�

%��ȡ����Phi������������,Mû����
[M,N]=size(Phi);

%�жϾ���(Phi'*Phi)�Ƿ����
if(rank(Phi'*Phi)~=N)
    error('���󲻿���');
end

y_p=Phi*((Phi'*Phi)\eye(N))*Phi'*y;
y_r=y-y_p;
end
