function [childp,childm,childf,F]=TSelection(p_chrom,m_chrom,f_chrom,fitness)
%�����Ƚ�ѡ����õĸ�����Ϊ�������뽻��أ�ֱ��ѡ��һ�ٸ����������������ظ�
    global ps SH N;
    pool_size=ps;%������ѡ��Ĳ���������صĴ�С
    tour=2;%���뾺���������Ŀ
    childp=zeros(ps,SH);
    childm=zeros(ps,SH);
    childf=zeros(ps,N);
    F=zeros(ps,3);
    for i=1:pool_size
        index1=ceil(rand*ps);
        index2=ceil(rand*ps);%���ѡ������λ��
        while index1==index2
            index2=ceil(rand*ps);
        end
        f1=fitness(index1,1:2);
        f2=fitness(index2,1:2);
        if(NDS(f1,f2)==1)
            childp(i,:)=p_chrom(index1,:);
            childm(i,:)=m_chrom(index1,:);
            childf(i,:)=f_chrom(index1,:);
            F(i,:)=fitness(index1,:);
        elseif(NDS(f1,f2)==2)
            childp(i,:)=p_chrom(index2,:);
            childm(i,:)=m_chrom(index2,:);
            childf(i,:)=f_chrom(index2,:);
            F(i,:)=fitness(index2,:);
        else
            if rand<=0.5
                childp(i,:)=p_chrom(index1,:);
                childm(i,:)=m_chrom(index1,:);
                childf(i,:)=f_chrom(index1,:);
                F(i,:)=fitness(index1,:);
            else
                childp(i,:)=p_chrom(index2,:);
                childm(i,:)=m_chrom(index2,:);
                childf(i,:)=f_chrom(index2,:);
                F(i,:)=fitness(index2,:);
            end
        end
    end
end