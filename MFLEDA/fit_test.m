function []=fit_test(p_chrom,m_chrom,f_chrom)
global N SH F;
%fit1:Cmax
%fit2�����ܺ�
%fit3���ؼ�����[0��F-1]
s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
    s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
end

P=cell(1,F);
for i=1:SH
    t1=s1(i);%��¼����ǰ���Ǹ�����
    t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
    P{f_chrom(t1)+1}=[P{f_chrom(t1)+1} p_chrom(i)];
end
FJ=cell(1,F);
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end

for i=1:F
    fitFJSP(P{i},m_chrom,FJ{i},i);
end



end

function fitFJSP(pro_m,mac_m,FJ,F_index)%���ô�ͳ��˫����뷽ʽ ������깤ʱ��makespan ��ʾ���һ����ɹ����ʱ��
global N H  TM time s_time;

Ep=4;%�ӹ����� 500kwÿs 0.5Mw
Ei=1;%�ȴ����� 40kwÿs
Es=2;
%     Ec=0.5;%���������ȹ��õĹ��� 60kwÿs

JOBN=length(FJ);
SH=length(pro_m);

finish={};%�����깤ʱ��
for i=1:JOBN%��ʼ�����ʱ�����
    JOBI=FJ(i);
    for j=1:H(JOBI)
        finish{i,j}=0;
    end
end

mt=cell(1,TM);
for i=1:TM%��ʼ������������ʱ������
    mt{i}=0;
end

s1=pro_m;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
    s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
end

for i=1:SH
    t1=s1(i);%��¼����ǰ���Ǹ�����
    t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
    mm(i)=mac_m(1,sum(H(1,1:t1-1))+t2);%��ȡ�ù���ôμӹ��Ļ���ѡ����Ϊ����������б�ʾ�ù����ڼ��μӹ���ѡ�Ļ�������һ�α�ʾһ������
end
%% �ж�

for i=1:SH
    if isempty(time{F_index,s1(i),s2(i),mm(i)})
        fprintf('\n%d\t%d\t%d\t%d\n%s\n',F_index,s1(i),s2(i),mm(i),'WRONG!');
    end
end
%% 



end