function new_p_chrom=GMEDA_p(p_chrom)
%��p_chrom��д GM-EDA�Ľ�ģ��������д
global ps N SH;
s_chrom=zeros(ps,SH);
s1=p_chrom(1,:);
p=zeros(1,N);
for i=1:SH
    p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
end

%%
%%%%%%��д
%2018-Solving the flexible job shop scheduling problem with sequence-dependent setup times
%Pseudocode 1(b):initialization for operation sequence vector
%p:������-->total operations per job
%p_chrom:tentative sequence vector
for pop_index=1:ps
    total_p=p;
    operation=1;
    job=1;
    while job<=N
        index=1;
        while index<=SH
            if p_chrom(pop_index,index)==job
                if total_p(job)>0
                    s_chrom(pop_index,index)=operation;
                    total_p(job)=total_p(job)-1;
                    operation=operation+1;
                end
            end
            index=index+1;
        end
        job=job+1;
    end
end

%%
%%%%%%GM-EDA��ģ������
%%%%%%����ѡ��ĸ��彨��ģ��%%%%%%
initialTheta=0.001;
upperTheta=10;
maxit=100;%����������
%1.���ݾ�����С���0
distances=zeros(ps);
for i=1:ps
    for j=i+1:ps
        perm1=s_chrom(i,:);
        perm2=s_chrom(j,:);
        [~,invperm2]=sort(perm2);
        composition1=perm1(invperm2);
        dist=sum(vVector(composition1));
        distances(i,j)=dist;
    end
end
distances=distances+distances';
sumdist=sum(distances);
[~,index]=min(sumdist);
ConsensusRanking=s_chrom(index,:);
%2.Newton-Raphson�㷨���j
[~,invCR]=sort(ConsensusRanking);
composition=s_chrom(:,invCR);
for i=1:ps
    vjs(i,:)=vVector(composition(i,:));
end
vjsmean=mean(vjs);
for j=1:SH-1
    Theta(j) = NewtonRaphson(initialTheta,vjsmean,upperTheta,maxit,@GKendallThetaFunction,@GKendallThetaDevFunction,{SH,j});
end
%3.���һ��������
j=1:SH-1;
Psi = (1 - exp((-1)*(SH-j+1).*Theta)) ./ (1 - exp((-1).*Theta));%��ʽ��7��
%4.��P(Vj(�Ҧ�0-1)=rj) ��ʽ��6��
for j=1:SH-1
    for r=0:SH-j
        upper=exp((-1)*r*Theta(j));
        lower=Psi(j);
        VProbs(j,r+1)=upper/lower;%ÿһ�к�Ϊ1��NumbVar-1�У���Ӧ���䵽NumbVar-1��λ�õĸ���
    end
end

%%
%%%%%%������������Ⱥ%%%%%%
NumbVar=SH;
NewPopSize=ps;
aux_v=zeros(1,NumbVar);
NewPop=zeros(NewPopSize,NumbVar);
randValues=rand(NewPopSize,NumbVar-1);
limit=NumbVar-1:-1:1;
for j=1:NewPopSize
    for i=1:NumbVar-1
        randVal=randValues(j,i);
        accumul=VProbs(i,1);
        index=0;
        while (index<limit(i))&&(accumul<randVal)
            accumul=accumul+VProbs(i,(index+1)+1);
            index=index+1;
        end
        aux_v(i)=index;
    end
    aux_v(NumbVar)=0;
    %��v���µ�����
    %function [permu] = GeneratePermuFromV(v,NumbVar)
    aux_n=1:NumbVar;
    for i=1:NumbVar-1
        val=aux_v(i);
        index=1;
        while (~(aux_n(index)~=-1&&val==0))
            %����1��aux_n(index)~=-1��val=0���˳�ѭ��
            %����2��aux_n(index)~=-1��val~=0��val-��index+
            %����3��aux_n(index)=-1��val=0��index+������index=1,2,3,...,ֱ��aux_n(index)~=-1����������1
            %����4��aux_n(index)=-1��val~=0��index+������index=1,2,3,...,��aux_n(index)~=-1��val-��index+
            if aux_n(index)~=-1
                val=val-1;
            end
            index=index+1;
        end
        permu(i)=index;%index�Ŀ���ȡֵ��[0,1,2,...,NumbVar-i]
        aux_n(index)=-1;
    end
    index=1;
    while (aux_n(index)==-1)
        index=index+1;
    end
    permu(NumbVar)=index;
    perm1=permu;
    perm2=ConsensusRanking;
    newpermutation=perm1(perm2);
    NewPop(j,1:NumbVar)=newpermutation;
end
p_chrom=NewPop;

%%
%%%%%%��д��ȥ
new_p_chrom=zeros(ps,SH);
for pop_index=1:ps
    total_p=p;
    new_chrom=p_chrom(pop_index,:);
    operation=1;
    for i=1:N
        while total_p(i)>0
            for j=1:SH                
                if new_chrom(j)==operation
                    new_p_chrom(pop_index,j)=i;
                    total_p(i)=total_p(i)-1;
                    operation=operation+1;
                    if total_p(i)<=0
                        break;
                    end
                end
            end
        end
    end
end

end

function newvector=vVector(vector)
%����1��29
for i=1:length(vector)-1
    newvector(i)=sum(vector(i+1:end)<vector(i));
end
end