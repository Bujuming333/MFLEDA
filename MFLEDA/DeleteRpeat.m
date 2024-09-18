function DeleteRpeat()
global AP AM AF AN;
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
            j=j-1;
            m=m-1;
        end
        j=j+1;
    end
end
end