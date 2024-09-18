function [TopPSRank]=FastNDS(fitness)
TopPSRank=[];
global ps;
L=size(fitness,1);
rank=1;
F(rank).f=[];     %��¼pareto�⼯�ȼ�Ϊrank���ĸ��弯��
individual=[];  %���ڴ�ű�ĳ������֧��ĸ��弯��
for i=1:L
    individual(i).n=0;  %n�Ǹ���i��֧��ĸ�������
    individual(i).p=[]; %p�Ǳ�����i֧��ĸ��弯��
end
for i=1:L
    for j=1:L
        dom_less=0;
        dom_equal=0;
        dom_more=0;
        for k=1:2  %�жϸ���i��j��֧���ϵ
            if (fitness(i,k)>fitness(j,k))
                dom_more = dom_more + 1;
            elseif (fitness(i,k)==fitness(j,k))
                dom_equal = dom_equal + 1;
            else
                dom_less = dom_less + 1;
            end
        end
        
        if dom_less == 0 && dom_equal ~= 2 % ˵��i��j֧�䣬��Ӧ��n��1
            individual(i).n = individual(i).n + 1;
        elseif dom_more == 0 && dom_equal ~= 2 % ˵��i֧��j,��j����i��֧��ϼ���
            individual(i).p = [individual(i).p j];
        end
    end
    if individual(i).n == 0 %����i��֧��ȼ�������ߣ����ڵ�ǰ���Ž⼯����Ӧ��Ⱦɫ����Я����������������Ϣ
        F(rank).f = [F(rank).f i];%�ȼ�Ϊ1�ķ�֧��⼯
    end
end
%����Ĵ�����Ϊ���ҳ��ȼ���ߵķ�֧��⼯
%����Ĵ�����Ϊ�˸�����������зּ�
while(~isempty(F(rank).f))
    Q=[];               %�����һ��front����
    fL=length(F(rank).f);
    for i=1:fL              %ѭ����ǰ֧��⼯�еĸ���
        if ~isempty(individual(F(rank).f(i)).p)
            pL=length(individual(F(rank).f(i)).p);
            for j=1:pL          %ѭ������i��֧��⼯�еĸ���
                k=individual(F(rank).f(i)).p(j);%��pf��i����֧���j��Ӧ����֧������1
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
%obj=fitness;%���и����Ŀ�꺯��ֵ��������ת���ɵ�ֵ�������㣬���������˵������ܽ�ӵ���Ⱥܴ����Կ���
%%FastNDS.m59-�޸�61-64�У�Ӧ�ð���������ǰ������fitness��
obj=[];
for front=1:(fL-1)
    obj=[obj;fitness(F(front).f,:)];
end
current_index=0;
for front=1:(fL-1)%F���һ��Ԫ��Ϊ��
    previous_index=current_index+1;%���ڱ��һ��ǰ��rank�ȼ��������ʼ��
    ffL=length(F(front).f);
    y=[];
    for i=1:ffL
        y(i,:)=obj(current_index+i,:);
    end
    
    current_index=current_index+i;
    crowd=[];
    crowd=zeros(ffL,2);%���ڴ������Ŀ��ֵ��ӵ����
    for i=1:2%%����Ŀ��ֵ
        [sort_based_on_objective,sort_index]=sort(y(:,i));%��һ��ǰ�����Ŀ��ֵ����������
        fmin=sort_based_on_objective(1);
        fmax=sort_based_on_objective(ffL);
        if ffL==1||fmax==fmin%������������Ĵ���
            for j=1:ffL
                crowd(j,i)=1;
            end
        else
            for j=1:ffL%����ÿ������ú���Ŀ��ֵ��ӵ����
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
    %ӵ����Խ��Խ�� ԽС˵���¶��� ӵ���ȴӴ�С�������� ͬʱ ��Ⱥ���������Ҳ��������
    tmp=[];sort_index=[];
    tmp=F(front).f;
    [~,sort_index]=sort(crowd(:,1),'descend');%ӵ���Ƚ�������
    %         fprintf('%d\r\n',front);
    if ffL<length(sort_index)
        fprintf('%d %d\r\n',ffL,length(sort_index));
    end
    for i=1:length(sort_index)
        F(front).f(i)=tmp(sort_index(i));
    end
end

%% ѡǰ���100������
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