function [newp,newm,newf]=LS2(p_chrom,m_chrom,f_chrom,fitness)
%随机机器选择
global N SH F H NM M;

s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
    s2(i)=p(s1(i));%记录加工过程中，工件的次数
end

P=cell(1,F);IP=cell(1,F);%P：每个工厂的加工工序；IP：每个工厂加工工序索引
for i=1:SH
    t1=s1(i);%记录到当前是那个工件
    t2=s2(i);%记录当前工件是加工到第几次
    P{f_chrom(t1)+1}=[P{f_chrom(t1)+1} p_chrom(i)];
    IP{f_chrom(t1)+1}=[IP{f_chrom(t1)+1} i];
end

FJ=cell(1,F);%FJ：每个工厂加工的工件号
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end

%关键工厂
critical_f=fitness(3)+1;

%随机选择一个工序
L=length(P{critical_f});%L：关键工厂加工的总工序数
Index=ceil(rand*L);
Job=s1(IP{critical_f}(Index));
Op=s2(IP{critical_f}(Index));
while NM{Job,Op}==1
    Index=ceil(rand*L);
    Job=s1(IP{critical_f}(Index));
    Op=s2(IP{critical_f}(Index));
end

%%%调整加工机器
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