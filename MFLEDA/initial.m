function [p_chrom,m_chrom,f_chrom] = initial()%该函数采用完全随机初始化的方法以及传统双层编码形式。
global  N F H SH NM ps M;%增加F 工厂数-->工厂码0~F-1

m_chrom=zeros(ps,SH);%机器码
p_chrom=zeros(ps,SH);%工序码
f_chrom=zeros(ps,N);%工厂码
chrom=zeros(1,SH);
FChrom=zeros(1,N);
%生成工序码
for i=1:N%将工序按顺序排开
    for j=1:H(i)
        k=sum(H(1,1:i))-H(i)+j;
        chrom(k)=i;
    end
    FChrom(i)=mod(i,F);%0表示分到工厂1,1表示分到工厂2,……F-1表示分到工厂F
end
%生成第一个工序码、工厂码
tmp=chrom;
p_chrom(1,:)=tmp(randperm(length(tmp)));%将工序码打乱
tmp2=FChrom;
f_chrom(1,:)=tmp2(randperm(length(tmp2)));%将工厂码打乱

for i=2:ps
    tmp=p_chrom(i-1,:);
    p_chrom(i,:)=tmp(randperm(length(tmp)));%再根据上一个生成后续的种群个体
    tmp2=f_chrom(i-1,:);
    f_chrom(i,:)=tmp2(randperm(length(tmp2)));%再根据上一个生成后续的种群个体
end
%工序码种群生成完毕

%生成机器码
for k=1:ps
    for i=1:N
        for j=1:H(i)
            t=ceil(rand*NM{i,j});
            t1=sum(H(1,1:i-1))+j;
            m_chrom(k,t1)=M{i,j,t};
        end
    end
end
%机器码种群生成完毕
end