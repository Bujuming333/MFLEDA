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
            j=j-1;%%Ϊʲô-1��AP,AM,AF,AN����һ�У���һ�б�Ϊ��ǰj����22��+1ǰҪ��һ
            m=m-1;
            if m<ps*2+1 %��֤����� 2*ps�����ѡ
                break;
            end
        end
        j=j+1;
    end
    if m<ps+1 %��֤�����ps�������ѡ
        break;
    end
end
end