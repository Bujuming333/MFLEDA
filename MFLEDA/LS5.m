function [newp,newm,newf]=LS5(p_chrom,m_chrom,f_chrom,fitness) %critical-block based neighborhood 4--N4
%N6��
%�ҵ��ؼ�·���Լ�����ͬһ�������������������Ϊһ���ٽ�飬
%����ͷ���β����м�飬
%���ѡһ�����뵽ͷ������ǰ�����ѡһ�����뵽β�������
%ͷ���β��ģ���ͷ���β��������ǰ�Ĺ�����뵽β�������
%��β��ͷ������ĺ���Ĺ�����뵽ͷ������ǰ��

global N SH F;

s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

newf=f_chrom;

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

for i=1:F
    if fitness(1,3)==i-1%�ؼ�����
        [CriticalPath,CriticalBlock,block]=FindCriticalPathDFJSP(P{i},m_chrom,FJ{i},i); %�ؼ�·�����ص�������Ⱦɫ���йؼ�������±�
    end
end

for i=1:block
    BL=length(CriticalBlock(i).B);
    if BL>1
        
        if i==1 %�׸��ٽ�죺���ѡһ��������뵽β�������
            Index1=ceil(rand*(BL-1));
            Index2=BL;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            for f_index=1:F     %����P{}
               if  fitness(1,3)==f_index-1
                   tmp=P{f_index}(Index1);
                   for j=Index1:Index2-1
                       P{f_index}(j)=P{f_index}(j+1);
                   end
                   P{f_index}(Index2)=tmp;
               end                
            end
        end
        
        if i==block %���һ���ٽ�죺���ѡһ��������뵽ͷ������ǰ
            Index1=1;
            Index2=ceil(rand*(BL-1))+1;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            for f_index=1:F     %����P{}
               if  fitness(1,3)==f_index-1
                   tmp=P{f_index}(Index2);
                   for j=Index2:-1:Index1+1
                       P{f_index}(j)=P{f_index}(j-1);
                   end
                   P{f_index}(Index1)=tmp;
               end                
            end
        end
        
        if i>1&&i<block&&BL>2 %�м��ٽ�죺���ѡһ��������뵽ͷ������ǰ�����ѡһ�����뵽β�������
            Index1=ceil(rand*(BL-2))+1; %�м����뵽β��
            Index2=BL;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            for f_index=1:F     %����P{}
               if  fitness(1,3)==f_index-1
                   tmp=P{f_index}(Index1);
                   for j=Index1:Index2-1
                       P{f_index}(j)=P{f_index}(j+1);
                   end
                   P{f_index}(Index2)=tmp;
               end                
            end
            
            Index1=1; %�м�����ͷ��֮ǰ
            Index2=ceil(rand*(BL-2))+1;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            for f_index=1:F     %����P{}
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