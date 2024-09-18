function [TopPSRank]=FastNDS(fitness)
TopPSRank=[];
global ps;
L=size(fitness,1);
rank=1;
F(rank).f=[];     %记录pareto解集等级为rank级的个体集合
individual=[];  %用于存放被某个个体支配的个体集合
for i=1:L
    individual(i).n=0;  %n是个体i被支配的个体数量
    individual(i).p=[]; %p是被个体i支配的个体集合
end
for i=1:L
    for j=1:L
        dom_less=0;
        dom_equal=0;
        dom_more=0;
        for k=1:2  %判断个体i和j的支配关系
            if (fitness(i,k)>fitness(j,k))
                dom_more = dom_more + 1;
            elseif (fitness(i,k)==fitness(j,k))
                dom_equal = dom_equal + 1;
            else
                dom_less = dom_less + 1;
            end
        end
        
        if dom_less == 0 && dom_equal ~= 2 % 说明i受j支配，相应的n加1
            individual(i).n = individual(i).n + 1;
        elseif dom_more == 0 && dom_equal ~= 2 % 说明i支配j,把j加入i的支配合集中
            individual(i).p = [individual(i).p j];
        end
    end
    if individual(i).n == 0 %个体i非支配等级排序最高，属于当前最优解集，相应的染色体中携带代表排序数的信息
        F(rank).f = [F(rank).f i];%等级为1的非支配解集
    end
end
%上面的代码是为了找出等级最高的非支配解集
%下面的代码是为了给其他个体进行分级
while(~isempty(F(rank).f))
    Q=[];               %存放下一个front集合
    fL=length(F(rank).f);
    for i=1:fL              %循环当前支配解集中的个体
        if ~isempty(individual(F(rank).f(i)).p)
            pL=length(individual(F(rank).f(i)).p);
            for j=1:pL          %循环个体i所支配解集中的个体
                k=individual(F(rank).f(i)).p(j);%被pf中i个体支配的j对应的受支配数减1
                individual(k).n=individual(k).n-1;
                if individual(k).n==0
                    Q=[Q k];
                end
            end
        end
    end
    rank=rank+1;
    F(rank).f=Q;
end

fL=length(F);
%obj=fitness;%所有个体的目标函数值经过处理转化成单值进行运算，就算结果相等说明距离很近拥挤度很大，所以可行
%%FastNDS.m59-修改61-64行：应该按照帕累托前沿排序fitness后
obj=[];
for front=1:(fL-1)
    obj=[obj;fitness(F(front).f,:)];
end
current_index=0;
for front=1:(fL-1)%F最后一个元素为空
    previous_index=current_index+1;%用于标记一个前沿rank等级个体的起始点
    ffL=length(F(front).f);
    y=[];
    for i=1:ffL
        y(i,:)=obj(current_index+i,:);
    end
    
    current_index=current_index+i;
    crowd=[];
    crowd=zeros(ffL,2);%用于存放两个目标值的拥挤度
    for i=1:2%%两个目标值
        [sort_based_on_objective,sort_index]=sort(y(:,i));%把一个前沿面的目标值做升序排序
        fmin=sort_based_on_objective(1);
        fmax=sort_based_on_objective(ffL);
        if ffL==1||fmax==fmin%对于特殊情况的处理
            for j=1:ffL
                crowd(j,i)=1;
            end
        else
            for j=1:ffL%计算每个个体该函数目标值的拥挤度
                if j==1
                    crowd(sort_index(j),i)=(sort_based_on_objective(2)-sort_based_on_objective(1))/(fmax-fmin);
                elseif j==ffL
                    crowd(sort_index(j),i)=(sort_based_on_objective(j)-sort_based_on_objective(j-1))/(fmax-fmin);
                else
                    crowd(sort_index(j),i)=(sort_based_on_objective(j+1)-sort_based_on_objective(j-1))/(fmax-fmin);
                end
            end
        end
        
    end
    crowd(:,1)=crowd(:,1)+crowd(:,2);
    %拥挤度越大越好 越小说明月堵塞 拥挤度从大到小降序排序 同时 种群个体的索引也跟着排序
    tmp=[];sort_index=[];
    tmp=F(front).f;
    [~,sort_index]=sort(crowd(:,1),'descend');%拥挤度降序排序
    %         fprintf('%d\r\n',front);
    if ffL<length(sort_index)
        fprintf('%d %d\r\n',ffL,length(sort_index));
    end
    for i=1:length(sort_index)
        F(front).f(i)=tmp(sort_index(i));
    end
end

%% 选前面的100个个体
count=0;
for front=1:(fL-1)
    ffL=length(F(front).f);
    for i=1:ffL
        TopPSRank=[TopPSRank,F(front).f(i)];
        count=count+1;
        if count==ps
            break;
        end
    end
    if count==ps
        break;
    end
end

end