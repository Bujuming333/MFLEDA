function [PF,PF1]=pareto1(obj)%根据两个适应函数 求出非支配解集
    global ps;
    PF=[];%AS存储种群中非支配解集的下标
    PF1=[];%AS存储种群中非支配解集的下标
    M=2;
    [obj_size,~]=size(obj);
    pn=zeros(1,obj_size);
    Nim=6;
    PF1=[];
    S=0;
    for i=1:obj_size
        for j=1:obj_size
            dom_less=0;
            dom_equal=0;
            dom_more=0;
            if (obj(i,1)>obj(j,1))  
                    dom_more = dom_more + 1;
            elseif (obj(i,1)==obj(j,1))
                    dom_equal = dom_equal + 1;
            else
                    dom_less = dom_less + 1;
            end
            
            if (obj(i,2)>obj(j,2))  
                    dom_more = dom_more + 1;
            elseif (obj(i,2)==obj(j,2))
                    dom_equal = dom_equal + 1;
            else
                    dom_less = dom_less + 1;
            end
            
            if dom_less == 0 && dom_equal ~= M % 说明i受j支配，相应的n加1
                pn(i) = pn(i)+ 1;
            end
        end
        
        if pn(i)== 0 %个体i非支配等级排序最高，属于当前最优解集，相应的染色体中携带代表排序数的信息
            PF=[PF i];
            S=S+1;
        end
    end
        


end