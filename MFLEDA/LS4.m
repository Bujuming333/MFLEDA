function [newp,newm,newf]=LS4(p_chrom,m_chrom,f_chrom,fitness)%one point insert, random insert neighborhood
%�ڹؼ�����insert,ÿ������ѡ��Ļ�������
%����������
global N H time SH F;%

s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
    s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
end

P=cell(1,F);IP=cell(1,F);
for i=1:SH
    t1=s1(i);%��¼����ǰ���Ǹ�����
    t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
    P{f_chrom(t1)+1}=[P{f_chrom(t1)+1} p_chrom(i)];
    IP{f_chrom(t1)+1}=[IP{f_chrom(t1)+1} i];
end

FJ=cell(1,F);
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end

%�ؼ�����
critical_f=fitness(3)+1;
L=length(P{critical_f});

IndexO1=ceil(rand*L);IndexO2=ceil(rand*L);
while IndexO1==IndexO2
    IndexO2=ceil(rand*L);
end
if IndexO1>IndexO2
    tmp=IndexO1;
    IndexO1=IndexO2;
    IndexO2=tmp;
end

%���ݼӹ�˳��IP{critical_f}��������֤ͬһ����

%%%��������
newp=p_chrom;
tmp=newp(IP{critical_f}(IndexO2));%ȡ��t2���ֵ�ֵ
%��t1��t2-1����ֵ�������
for i=IndexO2:-1:IndexO1+1
    newp(IP{critical_f}(i))=newp(IP{critical_f}(i-1));
end
newp(IP{critical_f}(IndexO1))=tmp;%��t2��ֵ���뵽t1��

%%%�����ӹ�����
newm=m_chrom;

%%%����s2
news2=s2;
tmp=news2(IP{critical_f}(IndexO2));%ȡ��t2���ֵ�ֵ
%��t1��t2-1����ֵ�������
for i=IndexO2:-1:IndexO1+1
    news2(IP{critical_f}(i))=news2(IP{critical_f}(i-1));
end
news2(IP{critical_f}(IndexO1))=tmp;%��t2��ֵ���뵽t1��

newf=f_chrom;
end