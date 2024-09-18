function [newp,newm,newf]=LS4(p_chrom,m_chrom,f_chrom,fitness)%one point insert, random insert neighborhood
%在关键工厂insert,每个工序选择的机器不变
%随机单点插入
global N H time SH F;%

s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
    s2(i)=p(s1(i));%记录加工过程中，工件的次数
end

P=cell(1,F);IP=cell(1,F);
for i=1:SH
    t1=s1(i);%记录到当前是那个工件
    t2=s2(i);%记录当前工件是加工到第几次
    P{f_chrom(t1)+1}=[P{f_chrom(t1)+1} p_chrom(i)];
    IP{f_chrom(t1)+1}=[IP{f_chrom(t1)+1} i];
end

FJ=cell(1,F);
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end

%关键工厂
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

%根据加工顺序IP{critical_f}调整，保证同一工厂

%%%调整工序
newp=p_chrom;
tmp=newp(IP{critical_f}(IndexO2));%取出t2部分的值
%将t1到t2-1部分值整体后移
for i=IndexO2:-1:IndexO1+1
    newp(IP{critical_f}(i))=newp(IP{critical_f}(i-1));
end
newp(IP{critical_f}(IndexO1))=tmp;%将t2的值插入到t1处

%%%调整加工机器
newm=m_chrom;

%%%调整s2
news2=s2;
tmp=news2(IP{critical_f}(IndexO2));%取出t2部分的值
%将t1到t2-1部分值整体后移
for i=IndexO2:-1:IndexO1+1
    news2(IP{critical_f}(i))=news2(IP{critical_f}(i-1));
end
news2(IP{critical_f}(IndexO1))=tmp;%将t2的值插入到t1处

newf=f_chrom;
end