function [f_chrom]=UMDA_f(Par_f)
global  N ps F;
f_chrom=zeros(ps,N);
%对N个工件：
%F个工厂，Par_m：[0,F-1]，分配到工厂1……工厂F
%1.求个数
I=zeros(F,N);
for i=1:N
    for j=1:ps
        for k=1:F
            if Par_f(j,i)==k-1        %工厂k
                I(k,i)=I(k,i)+1;
            end
        end
    end
end
%2.求概率
prob1=I/ps;
prob1=cumsum(prob1);
%3.生成新解
for i=1:ps
    for j=1:N
        r=rand;
        for k=1:F            
            if r<prob1(k,j)
                f_chrom(i,j)=k-1;   %工厂k
                break;
            end
        end
    end
end