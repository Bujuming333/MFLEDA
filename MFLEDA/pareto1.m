function [PF,PF1]=pareto1(obj)%����������Ӧ���� �����֧��⼯
    global ps;
    PF=[];%AS�洢��Ⱥ�з�֧��⼯���±�
    PF1=[];%AS�洢��Ⱥ�з�֧��⼯���±�
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
            
            if dom_less == 0 && dom_equal ~= M % ˵��i��j֧�䣬��Ӧ��n��1
                pn(i) = pn(i)+ 1;
            end
        end
        
        if pn(i)== 0 %����i��֧��ȼ�������ߣ����ڵ�ǰ���Ž⼯����Ӧ��Ⱦɫ����Я����������������Ϣ
            PF=[PF i];
            S=S+1;
        end
    end
        


end