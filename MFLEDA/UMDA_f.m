function [f_chrom]=UMDA_f(Par_f)
global  N ps F;
f_chrom=zeros(ps,N);
%��N��������
%F��������Par_m��[0,F-1]�����䵽����1��������F
%1.�����
I=zeros(F,N);
for i=1:N
    for j=1:ps
        for k=1:F
            if Par_f(j,i)==k-1        %����k
                I(k,i)=I(k,i)+1;
            end
        end
    end
end
%2.�����
prob1=I/ps;
prob1=cumsum(prob1);
%3.�����½�
for i=1:ps
    for j=1:N
        r=rand;
        for k=1:F            
            if r<prob1(k,j)
                f_chrom(i,j)=k-1;   %����k
                break;
            end
        end
    end
end