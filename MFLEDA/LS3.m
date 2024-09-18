function [newp,newm,newf]=LS3(p_chrom,m_chrom,f_chrom,fitness)%swap Two-point exchange neigborhood
%在关键工厂随机交换,每个工序选择的机器不变
global N SH F;

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
L=length(P{critical_f});%L：关键工厂加工的总工序数

IndexO1=ceil(rand*L);IndexO2=ceil(rand*L);
while IndexO1==IndexO2
    IndexO2=ceil(rand*L);
end

%%%调整加工工序
newp=p_chrom;
tmp=newp(IP{critical_f}(IndexO1));
newp(IP{critical_f}(IndexO1))=newp(IP{critical_f}(IndexO2));%交换两点的染色体值
newp(IP{critical_f}(IndexO2))=tmp;

newm=m_chrom;

newf=f_chrom;
end




% tmp=P{critical_f}(IndexO1);
% P{critical_f}(IndexO1)=P{critical_f}(IndexO2);%交换两点的染色体值
% P{critical_f}(IndexO2)=tmp;
% 
% newp=p_chrom;
% for i=1:L
%     newp(1,IP{critical_f}(i))=P{critical_f}(i);
% end