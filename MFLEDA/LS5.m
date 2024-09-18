function [newp,newm,newf]=LS5(p_chrom,m_chrom,f_chrom,fitness) %critical-block based neighborhood 4--N4
%N6：
%找到关键路径以及将在同一机器的连续操作工序归为一个临界块，
%不是头块和尾块的中间块，
%随机选一个插入到头部工序前，随机选一个插入到尾部工序后，
%头块和尾块的，将头块的尾部工序以前的工序插入到尾部工序后，
%将尾块头部工序的后面的工序插入到头部工序前。

global N SH F;

s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

newf=f_chrom;

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

for i=1:F
    if fitness(1,3)==i-1%关键工厂
        [CriticalPath,CriticalBlock,block]=FindCriticalPathDFJSP(P{i},m_chrom,FJ{i},i); %关键路径返回的是在子染色体中关键工序的下标
    end
end

for i=1:block
    BL=length(CriticalBlock(i).B);
    if BL>1
        
        if i==1 %首个临界快：随机选一个工序插入到尾部工序后
            Index1=ceil(rand*(BL-1));
            Index2=BL;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            for f_index=1:F     %更新P{}
               if  fitness(1,3)==f_index-1
                   tmp=P{f_index}(Index1);
                   for j=Index1:Index2-1
                       P{f_index}(j)=P{f_index}(j+1);
                   end
                   P{f_index}(Index2)=tmp;
               end                
            end
        end
        
        if i==block %最后一个临界快：随机选一个工序插入到头部工序前
            Index1=1;
            Index2=ceil(rand*(BL-1))+1;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            for f_index=1:F     %更新P{}
               if  fitness(1,3)==f_index-1
                   tmp=P{f_index}(Index2);
                   for j=Index2:-1:Index1+1
                       P{f_index}(j)=P{f_index}(j-1);
                   end
                   P{f_index}(Index1)=tmp;
               end                
            end
        end
        
        if i>1&&i<block&&BL>2 %中间临界快：随机选一个工序插入到头部工序前和随机选一个插入到尾部工序后
            Index1=ceil(rand*(BL-2))+1; %中间块插入到尾块
            Index2=BL;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            for f_index=1:F     %更新P{}
               if  fitness(1,3)==f_index-1
                   tmp=P{f_index}(Index1);
                   for j=Index1:Index2-1
                       P{f_index}(j)=P{f_index}(j+1);
                   end
                   P{f_index}(Index2)=tmp;
               end                
            end
            
            Index1=1; %中间块插入头部之前
            Index2=ceil(rand*(BL-2))+1;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            for f_index=1:F     %更新P{}
                if  fitness(1,3)==f_index-1
                    tmp=P{f_index}(Index2);
                    for j=Index2:-1:Index1+1
                        P{f_index}(j)=P{f_index}(j-1);
                    end
                    P{f_index}(Index1)=tmp;
                end
            end
        end
    end
end

newm=m_chrom;
newp=zeros(1,SH);
for f_index=1:F
    L=length(IP{f_index});
    for i=1:L
        newp(1,IP{f_index}(i))=P{f_index}(i);
    end
end

end