function [newp,newm,newf]=LS1(p_chrom,m_chrom,f_chrom,fitness)
%随机工厂分配
global N SH F;

FJ=cell(1,F);%FJ：每个工厂加工的工件号
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end

%关键工厂
critical_f=fitness(3)+1;
L=length(FJ{critical_f});
Index=ceil(rand*L);
Job=FJ{critical_f}(Index);
IndexF=ceil(rand*F);
while IndexF==critical_f
    IndexF=ceil(rand*F);
end
newf=f_chrom;
newf(Job)=IndexF-1;

newp=p_chrom;

newm=m_chrom;


end