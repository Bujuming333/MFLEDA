function new_p_chrom=GMEDA_p(p_chrom)
%对p_chrom改写 GM-EDA的建模、采样改写
global ps N SH;
s_chrom=zeros(ps,SH);
s1=p_chrom(1,:);
p=zeros(1,N);
for i=1:SH
    p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
end

%%
%%%%%%改写
%2018-Solving the flexible job shop scheduling problem with sequence-dependent setup times
%Pseudocode 1(b):initialization for operation sequence vector
%p:操作数-->total operations per job
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
%%%%%%GM-EDA建模、采样
%%%%%%根据选择的个体建立模型%%%%%%
initialTheta=0.001;
upperTheta=10;
maxit=100;%最大迭代次数
%1.根据距离最小求σ0
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
%2.Newton-Raphson算法求θj
[~,invCR]=sort(ConsensusRanking);
composition=s_chrom(:,invCR);
for i=1:ps
    vjs(i,:)=vVector(composition(i,:));
end
vjsmean=mean(vjs);
for j=1:SH-1
    Theta(j) = NewtonRaphson(initialTheta,vjsmean,upperTheta,maxit,@GKendallThetaFunction,@GKendallThetaDevFunction,{SH,j});
end
%3.求归一化常数ψ
j=1:SH-1;
Psi = (1 - exp((-1)*(SH-j+1).*Theta)) ./ (1 - exp((-1).*Theta));%公式（7）
%4.求P(Vj(σσ0-1)=rj) 公式（6）
for j=1:SH-1
    for r=0:SH-j
        upper=exp((-1)*r*Theta(j));
        lower=Psi(j);
        VProbs(j,r+1)=upper/lower;%每一行和为1，NumbVar-1行，对应分配到NumbVar-1个位置的概率
    end
end

%%
%%%%%%采样生成新种群%%%%%%
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
    %由v求新的排列
    %function [permu] = GeneratePermuFromV(v,NumbVar)
    aux_n=1:NumbVar;
    for i=1:NumbVar-1
        val=aux_v(i);
        index=1;
        while (~(aux_n(index)~=-1&&val==0))
            %情形1：aux_n(index)~=-1，val=0：退出循环
            %情形2：aux_n(index)~=-1，val~=0：val-，index+
            %情形3：aux_n(index)=-1，val=0：index+，遍历index=1,2,3,...,直到aux_n(index)~=-1，满足情形1
            %情形4：aux_n(index)=-1，val~=0：index+，遍历index=1,2,3,...,当aux_n(index)~=-1，val-，index+
            if aux_n(index)~=-1
                val=val-1;
            end
            index=index+1;
        end
        permu(i)=index;%index的可能取值：[0,1,2,...,NumbVar-i]
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
%%%%%%改写回去
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
%返回1×29
for i=1:length(vector)-1
    newvector(i)=sum(vector(i+1:end)<vector(i));
end
end