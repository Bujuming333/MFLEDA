function [newp,newm,newf]=LS3(p_chrom,m_chrom,f_chrom,fitness)%swap Two-point exchange neigborhood
%�ڹؼ������������,ÿ������ѡ��Ļ�������
global N SH F;

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
L=length(P{critical_f});%L���ؼ������ӹ����ܹ�����

IndexO1=ceil(rand*L);IndexO2=ceil(rand*L);
while IndexO1==IndexO2
    IndexO2=ceil(rand*L);
end

%%%�����ӹ�����
newp=p_chrom;
tmp=newp(IP{critical_f}(IndexO1));
newp(IP{critical_f}(IndexO1))=newp(IP{critical_f}(IndexO2));%���������Ⱦɫ��ֵ
newp(IP{critical_f}(IndexO2))=tmp;

newm=m_chrom;

newf=f_chrom;
end




% tmp=P{critical_f}(IndexO1);
% P{critical_f}(IndexO1)=P{critical_f}(IndexO2);%���������Ⱦɫ��ֵ
% P{critical_f}(IndexO2)=tmp;
% 
% newp=p_chrom;
% for i=1:L
%     newp(1,IP{critical_f}(i))=P{critical_f}(i);
% end