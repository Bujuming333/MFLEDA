function [newp,newm,newf]=LS2(p_chrom,m_chrom,f_chrom,fitness)
%�������ѡ��
global N SH F H NM M;

s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
    s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
end

P=cell(1,F);IP=cell(1,F);%P��ÿ�������ļӹ�����IP��ÿ�������ӹ���������
for i=1:SH
    t1=s1(i);%��¼����ǰ���Ǹ�����
    t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
    P{f_chrom(t1)+1}=[P{f_chrom(t1)+1} p_chrom(i)];
    IP{f_chrom(t1)+1}=[IP{f_chrom(t1)+1} i];
end

FJ=cell(1,F);%FJ��ÿ�������ӹ��Ĺ�����
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end

%�ؼ�����
critical_f=fitness(3)+1;

%���ѡ��һ������
L=length(P{critical_f});%L���ؼ������ӹ����ܹ�����
Index=ceil(rand*L);
Job=s1(IP{critical_f}(Index));
Op=s2(IP{critical_f}(Index));
while NM{Job,Op}==1
    Index=ceil(rand*L);
    Job=s1(IP{critical_f}(Index));
    Op=s2(IP{critical_f}(Index));
end

%%%�����ӹ�����
position=sum(H(1,1:Job-1))+Op;
IndexM=ceil(rand*NM{Job,Op});
while M{Job,Op,IndexM}==m_chrom(position)
    IndexM=ceil(rand*NM{Job,Op});
end
newm=m_chrom;
newm(position)=M{Job,Op,IndexM};

newp=p_chrom;

newf=f_chrom;
end