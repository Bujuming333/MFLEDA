function [AP,AM,AN,AF]=DeleteRpeatQF(AP,AM,AN,AF)
global ps;
[m,~]=size(AF);
for i=1:m
    if i>m
        break;
    end
    F=AF(i,:);
    j=i+1;
    while j<m
        if AF(j,1)==F(1,1)&&AF(j,2)==F(1,2)
            AP(j,:)=[];
            AM(j,:)=[];
            AF(j,:)=[];
            AN(j,:)=[];
            j=j-1;%%为什么-1：AP,AM,AF,AN少了一行，下一行变为当前j，第22行+1前要减一
            m=m-1;
            if m<ps*2+1 %保证最后有 2*ps个体参选
                break;
            end
        end
        j=j+1;
    end
    if m<ps+1 %保证最后有ps个个体参选
        break;
    end
end
end