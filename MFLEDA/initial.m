function [p_chrom,m_chrom,f_chrom] = initial()%�ú���������ȫ�����ʼ���ķ����Լ���ͳ˫�������ʽ��
global  N F H SH NM ps M;%����F ������-->������0~F-1

m_chrom=zeros(ps,SH);%������
p_chrom=zeros(ps,SH);%������
f_chrom=zeros(ps,N);%������
chrom=zeros(1,SH);
FChrom=zeros(1,N);
%���ɹ�����
for i=1:N%������˳���ſ�
    for j=1:H(i)
        k=sum(H(1,1:i))-H(i)+j;
        chrom(k)=i;
    end
    FChrom(i)=mod(i,F);%0��ʾ�ֵ�����1,1��ʾ�ֵ�����2,����F-1��ʾ�ֵ�����F
end
%���ɵ�һ�������롢������
tmp=chrom;
p_chrom(1,:)=tmp(randperm(length(tmp)));%�����������
tmp2=FChrom;
f_chrom(1,:)=tmp2(randperm(length(tmp2)));%�����������

for i=2:ps
    tmp=p_chrom(i-1,:);
    p_chrom(i,:)=tmp(randperm(length(tmp)));%�ٸ�����һ�����ɺ�������Ⱥ����
    tmp2=f_chrom(i-1,:);
    f_chrom(i,:)=tmp2(randperm(length(tmp2)));%�ٸ�����һ�����ɺ�������Ⱥ����
end
%��������Ⱥ�������

%���ɻ�����
for k=1:ps
    for i=1:N
        for j=1:H(i)
            t=ceil(rand*NM{i,j});
            t1=sum(H(1,1:i-1))+j;
            m_chrom(k,t1)=M{i,j,t};
        end
    end
end
%��������Ⱥ�������
end