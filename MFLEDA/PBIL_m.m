function [m_chrom,A]=PBIL_m(Par_m,A,alpha)
global  N H SH NM ps M;
m_chrom=zeros(ps,SH);%机器码100*55
SI=size(Par_m,1);
A0=A;
%求概率矩阵并返回
I=cell(1,N);
prob=I;
for n=1:N
    I{n}=zeros(H(n),5);
    for i=1:H(n)
        for j=1:SI
            for b=1:NM{n,i}
                if Par_m(j,sum(H(1,1:n-1))+i)==M{n,i,b}
                    I{n}(i,b)=I{n}(i,b)+1;
                end
            end
        end
    end
    A{n}=(1-alpha)*A{n}+alpha/SI*I{n};
    prob{n}=cumsum(A{n},2);
end
%轮盘赌选择加工机器
for k=1:ps
    for i=1:N
        for j=1:H(i)
            rand_r=rand;
            for x=1:NM{i,j}
                if rand_r<prob{i}(j,x)
                    t=x;
                    t1=sum(H(1,1:i-1))+j;
                    m_chrom(k,t1)=M{i,j,t};
                    break;
                end
            end
        end
    end
end
end