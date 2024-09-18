function [p_chrom,m_chrom,f_chrom,fitness]=EnergySaveDFJSP(p_chrom,m_chrom,f_chrom,fitness) %ͨ�����ٵȴ�ʱ����������������
global N SH AP AM AN AF F;

s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
    s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
end

P=cell(1,F);IP=cell(1,F);
for i=1:SH
    t1=s1(i);%��¼����ǰ���Ǹ�����
    t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
    P{f_chrom(t1)+1}=[P{f_chrom(t1)+1} p_chrom(i)];
    IP{f_chrom(t1)+1}=[IP{f_chrom(t1)+1} i];
end

FJ=cell(1,F);
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end

% drawFJSP_SDST_M5(P{1},m_chrom,FJ{1},1);
% drawFJSP_SDST_M5(P{2},m_chrom,FJ{2},2);
for i=1:F
    P{i}=SAS2AS(P{i},m_chrom,FJ{i},i);%��������
end
% drawFJSP_SDST_M5(P{1},m_chrom,FJ{1},1);
% drawFJSP_SDST_M5(P{2},m_chrom,FJ{2},2);
% for i=1:F
%     P{i}=AS2FAS(P{i},m_chrom,FJ{i},i);%ȫ��������
% end
% drawFJSP_SDST_M5(P{1},m_chrom,FJ{1},1);
% drawFJSP_SDST_M5(P{2},m_chrom,FJ{2},2);

new_f=f_chrom;
new_m=m_chrom;

new_p=zeros(1,SH);
for i=1:F
    L=length(IP{i});
    for j=1:L
        new_p(IP{i}(j))=P{i}(j);
    end
end

[newfit(1,1),newfit(1,2),newfit(1,3)]=fitDFJSP(new_p,new_m,new_f);
% [newfit2(1,1),newfit2(1,2),newfit2(1,3)]=fitDFJSP(new_p2,new_m,new_f);
if NDS(newfit,fitness)==1 %�½�֧��ɽ�
    p_chrom=new_p;
    m_chrom=new_m;
    f_chrom=new_f;
    fitness=newfit;
    AP=[AP;new_p];
    AM=[AM;new_m];
    AF=[AF;newfit];
    AN=[AN;new_f];
elseif NDS(newfit,fitness)==0 %�½⻥��֧��ɽ�
    AP=[AP;new_p];
    AM=[AM;new_m];
    AF=[AF;newfit];
    AN=[AN;new_f];
end

end

function [newp]=SAS2AS(p_chrom,m_chrom,FJ,f_index)%��������
global N H TM time s_time;

e=0;
JOBN=length(FJ);
SH=length(p_chrom);
finish={};%�����깤ʱ��
start={};%����ʼʱ��
for i=1:JOBN%��ʼ�����ʱ�����
    JOBI=FJ(i);
    for j=1:H(JOBI)
        finish{JOBI,j}=e;
        start{JOBI,j}=e;
    end
end

mt=cell(1,TM);
for i=1:TM%��ʼ������������ʱ������
    Machine(i).Op=[];
    Machine(i).GapT=[];
    Machine(i).MFT=[];
    mt{i}=e;
end
s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
    s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
end

for i=1:SH
    t1=s1(i);%��¼����ǰ���Ǹ�����
    t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
    mm(i)=m_chrom(1,sum(H(1,1:t1-1))+t2);%��ȡ�ù���ôμӹ��Ļ���ѡ����Ϊ����������б�ʾ�ù����ڼ��μӹ���ѡ�Ļ�������һ�α�ʾһ������
end
%��ʼ����
for i=1:SH
    if(s2(i)==1)
        ON=length(Machine(mm(i)).Op);%�û�����Ŀǰ�Ĺ�����
        if ON>0%ĳ�����ĵ�һ�����򣬵���ǰ����֮ǰ�мӹ�
            %��Ҫ�޸�t=time{f_index,s1(i),s2(i),mm(i)};�����Ǽӹ�ʱ��+׼��ʱ������бȽ�
            %t:1*ON
            t=zeros(1,ON);
            t(1)=time{f_index,s1(i),s2(i),mm(i)};%��һ������          
            if ON>1
                for ON_index=2:ON
                    %��ǰ�ӹ��Ĺ���s1(Machine(mm(i)).Op(ON_index-1))
                    t(ON_index)=time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(Machine(mm(i)).Op(ON_index-1)),s1(i));
                end
            end
            %127��
            Index1=0;
            for j=1:ON%��ǰ�����ҿռ�
                if Machine(mm(i)).GapT(j)-t(j)>0
                    Index1=j; %��Index����ǰ���������Ĺ���
                    break;
                end
            end
            if Index1~=0%�ɲ��룬���º������еļӹ��Ĺ�����Ϣ
                Index1=Machine(mm(i)).Op(Index1);% ��s1�����а�s1(i)���뵽s1(Index1)ǰ��,Index1<i
                
                %����s1
                tmp=s1(i);
                for j=i:-1:Index1+1
                    s1(j)=s1(j-1);
                end
                s1(Index1)=tmp;
                
                %����s2
                tmp=s2(i);
                for j=i:-1:Index1+1
                    s2(j)=s2(j-1);
                end
                s2(Index1)=tmp;
                
                %���¼ӹ�����mm��ͬ�ϣ�
                tmp=mm(i);%��ǰ����
                for j=i:-1:Index1+1
                    mm(j)=mm(j-1);
                end
                mm(Index1)=tmp;
                
                %���µ�ǰ�����ļӹ�����
                for j=1:ON
                    if Machine(mm(Index1)).Op(j)>=Index1
                        Machine(mm(Index1)).Op(j)=Machine(mm(Index1)).Op(j)+1;%����ǰ����빤������������Ҫ���������һ��
                    end
                end
                
                %�����������л����ļӹ�����
                for k=1:TM
                    if k~=mm(Index1)
                        ON2=length(Machine(k).Op);%�˿̻���k�ӹ����ܴ���
                        for h=1:ON2
                            if Machine(k).Op(h)>Index1&&Machine(k).Op(h)<i
                                Machine(k).Op(h)=Machine(k).Op(h)+1;
                            end
                        end
                    end
                end
                
                %���µ�ǰ�����ļӹ�������������Ĺ���
                Machine(mm(Index1)).Op=[Machine(mm(Index1)).Op,Index1];
                tmp2=Machine(mm(Index1)).Op;
                Machine(mm(Index1)).Op=sort(tmp2,'ascend');
                
                %���´����빤��Ŀ�ʼ�ӹ�ʱ��
                IIndex=find(Machine(mm(Index1)).Op==Index1);
                if IIndex==1%��Ϊ��ǰ�����ĵ�һ���ӹ�����ʼʱ��Ϊ0
                    start{s1(Index1),s2(Index1)}=0;
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t(1);
                else%���¿�ʼʱ�䣨��ʼ׼������max{�ù���ǰһ������Ľ���ʱ�䣬��ǰ����ǰһ������Ľ���ʱ��}
                    LastOp=Machine(mm(Index1)).Op(IIndex-1);%�����λ��ǰһ���ӹ�����
                    start{s1(Index1),s2(Index1)}=max(0,finish{s1(LastOp),s2(LastOp)});
                    %���´����빤��Ľ����ӹ�ʱ�䣻����׼��ʱ��
                    %��һ��������s1(LastOp)
                    t_newset=time{f_index,s1(Index1),s2(Index1),mm(Index1)}+s_time(s1(LastOp),s1(Index1));%׼��ʱ��+�ӹ�ʱ��
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t_newset;
                end
                
                %��¼���빤��
                insert_index=Index1;
                                
                %���µ�ǰ�����ļӹ���Ϣ
                ON=ON+1;
                for j=1:ON
                    Index1=Machine(mm(Index1)).Op(j);
                    if j==1
                        Machine(mm(Index1)).GapT(j)=0;
                    else%���ǵ�һ���ӹ�����ʱ
                        %�����ڲ��빤���ĵ�һ�����򣺿�ʼʱ�䲻�䣬׼��ʱ��ı�-->����ʱ��ı䣬ƽ��finish�Ĳ�ֵ
                        %�����ڲ��빤������������,���¹���Ŀ�ʼʱ�䡢���ʱ�䣺ƽ����ͬ�Ĳ�ֵ
                        %��������׼��ʱ�䣬���뵽����λ�ã�����Ĺ���ʼʱ�䲻��
                        LastOp=Machine(mm(Index1)).Op(j-1);%��һ������
                        if LastOp==insert_index%�����ڲ��빤���ĵ�һ������
                            new_finish=start{s1(Index1),s2(Index1)}+time{f_index,s1(Index1),s2(Index1),mm(Index1)}+s_time(s1(LastOp),s1(Index1));
                            shift_value=new_finish-finish{s1(Index1),s2(Index1)};
                            finish{s1(Index1),s2(Index1)}=new_finish;
                        end
                        if LastOp>insert_index%�����ڲ��빤�����������򣺿�ʼʱ�䡢���ʱ�䣺ƽ����ͬ�Ĳ�ֵ
                            start{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+shift_value;
                            finish{s1(Index1),s2(Index1)}=finish{s1(Index1),s2(Index1)}+shift_value;
                        end
                        %���¿���ʱ��
                        Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)}-finish{s1(LastOp),s2(LastOp)};
                        if Machine(mm(Index1)).GapT(j)<0
                            Machine(mm(Index1)).GapT(j)
                        end
                    end
                    %����ÿ��λ�õļӹ����ʱ��
                    Machine(mm(Index1)).MFT(j)=finish{s1(Index1),s2(Index1)};
                end
                
                %���»����깤ʱ��
                mt{mm(Index1)}=Machine(mm(Index1)).MFT(ON);
            else %��194�����index1==0˵��û�п�λ��Ҫ����ʵʵ��ȥ�ں���
                %��һ���ӹ�������s1(Machine(mm(i)).Op(end))
                start{s1(i),s2(i)}=Machine(mm(i)).MFT(ON);%Ҳ�ǿ�ʼ׼����ʱ��
                %�����깤ʱ��+׼��ʱ��
                mt{mm(i)}=start{s1(i),s2(i)}+time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(Machine(mm(i)).Op(end)),s1(i));
                finish{s1(i),s2(i)}=mt{mm(i)};
                Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
                Machine(mm(i)).GapT=[Machine(mm(i)).GapT,0];%�޿���
                Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
            end
            
        else%��203��ĳ�����ĵ�һ�������ҵ�ǰ������һ�μӹ������¹����Ŀ�ʼ������ʱ�䣬�����ӹ��Ĺ������ʱ��
            mt{mm(i)}=time{f_index,s1(i),s2(i),mm(i)};
            start{s1(i),s2(i)}=0;
            finish{s1(i),s2(i)}=mt{mm(i)};
            Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
            Machine(mm(i)).GapT=[Machine(mm(i)).GapT,0];%����ʱ��
            Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];%���ʱ��
        end
        
    else
        %213����ʼ����ͬ�Ļ�����Ѱ�Һ��ʵĿ�λ���룬�������ļӹ�ʱ��С�ڿ���ʱ������Բ���
        ON=length(Machine(mm(i)).Op);%�û�����Ŀǰ�Ĺ�����
        if ON>0
            %��Ҫ�޸�t=time{f_index,s1(i),s2(i),mm(i)};�����Ǽӹ�ʱ��+׼��ʱ������бȽ�
            %t:1*ON
            t=zeros(1,ON);
            t(1)=time{f_index,s1(i),s2(i),mm(i)};%��һ�����У���Ϊ��һ���ӹ�������׼��ʱ��          
            if ON>1
                for ON_index=2:ON
                    %��ǰ�ӹ��Ĺ���Machine(mm(i)).Op(ON_index-1)
                    t(ON_index)=time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(Machine(mm(i)).Op(ON_index-1)),s1(i));
                end
            end
                       
            Index1=0;
            for j=1:ON
                if Machine(mm(i)).GapT(j)>t(j)
                    if ON==1 || j==1
                        tmp=finish{s1(i),s2(i)-1}-0;
                    else
                        tmp=finish{s1(i),s2(i)-1}-Machine(mm(i)).MFT(j-1);
                    end
                    if Machine(mm(i)).GapT(j)-t(j)-tmp>0
                        Index1=j; %��Index����ǰ���������Ĺ���
                        break;
                    end
                end
            end
            if Index1~=0
                %232������ĳ�����ĵ�һ�����򣬿ɲ���
                Index1=Machine(mm(i)).Op(Index1);% ��s1�����а�s1(i)���뵽s1(Index1)ǰ��,Index1<i
                
                %����s1
                tmp=s1(i);
                for j=i:-1:Index1+1
                    s1(j)=s1(j-1);
                end
                s1(Index1)=tmp;
                
                %����s2
                tmp=s2(i);
                for j=i:-1:Index1+1
                    s2(j)=s2(j-1);
                end
                s2(Index1)=tmp;
                
                %���¼ӹ�����mm��ͬ�ϣ�
                tmp=mm(i);
                for j=i:-1:Index1+1
                    mm(j)=mm(j-1);
                end
                mm(Index1)=tmp;
                
                %���µ�ǰ�����ļӹ�����
                for j=1:ON
                    if Machine(mm(Index1)).Op(j)>=Index1
                        Machine(mm(Index1)).Op(j)=Machine(mm(Index1)).Op(j)+1;%����ǰ����빤������������Ҫ���������һ��
                    end
                end
                
                %�����������л����ļӹ�����
                for k=1:TM
                    if k~=mm(Index1)
                        ON2=length(Machine(k).Op);
                        for h=1:ON2
                            if Machine(k).Op(h)>Index1&&Machine(k).Op(h)<i
                                Machine(k).Op(h)=Machine(k).Op(h)+1;
                            end
                        end
                    end
                end
                
                %���µ�ǰ�����ļӹ�������������Ĺ���
                Machine(mm(Index1)).Op=[Machine(mm(Index1)).Op,Index1];
                tmp2=Machine(mm(Index1)).Op;
                Machine(mm(Index1)).Op=sort(tmp2,'ascend');
                
                %���´����빤��Ŀ�ʼ�ӹ�ʱ��
                IIndex=find(Machine(mm(Index1)).Op==Index1);
                if IIndex==1%��Ϊ��ǰ�����ĵ�һ���ӹ�����ʼʱ��Ϊ0
                    start{s1(Index1),s2(Index1)}=max(0,finish{s1(Index1),s2(Index1)-1});
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t(1);
                else%���¿�ʼʱ�䣨��ʼ׼������max{�ù���ǰһ������Ľ���ʱ�䣬��ǰ����ǰһ������Ľ���ʱ��}
                    LastOp=Machine(mm(Index1)).Op(IIndex-1);%�����λ��ǰһ���ӹ�����
                    start{s1(Index1),s2(Index1)}=max(finish{s1(Index1),s2(Index1)-1},finish{s1(LastOp),s2(LastOp)});
                    %���´����빤��Ľ����ӹ�ʱ�䣻����׼��ʱ��
                    %��һ��������s1(LastOp)
                    t_newset=time{f_index,s1(Index1),s2(Index1),mm(Index1)}+s_time(s1(LastOp),s1(Index1));%׼��ʱ��+�ӹ�ʱ��
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t_newset;
                end
                
                
                %��¼���빤��
                insert_index=Index1;
                
                %���µ�ǰ�����ļӹ���Ϣ
                ON=ON+1;
                shift_value=0;
                for j=1:ON                    
                    Index1=Machine(mm(Index1)).Op(j);
                    if j==1
                        Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)};
                    else%���ǵ�һ���ӹ�����ʱ
                        %�����ڲ��빤���ĵ�һ�����򣺿�ʼʱ�䲻�䣬׼��ʱ��ı�-->����ʱ��ı䣬ƽ��finish�Ĳ�ֵ
                        %�����ڲ��빤������������,���¹���Ŀ�ʼʱ�䡢���ʱ�䣺ƽ����ͬ�Ĳ�ֵ
                        %��������׼��ʱ�䣬���뵽����λ�ã�����Ĺ���ʼʱ�䲻��
                        LastOp=Machine(mm(Index1)).Op(j-1);%��һ������
                        if LastOp==insert_index%�����ڲ��빤���ĵ�һ������
                            new_finish=start{s1(Index1),s2(Index1)}+time{f_index,s1(Index1),s2(Index1),mm(Index1)}+s_time(s1(LastOp),s1(Index1));
                            shift_value=new_finish-finish{s1(Index1),s2(Index1)};
                            finish{s1(Index1),s2(Index1)}=new_finish;
                        end
                        if LastOp>insert_index%�����ڲ��빤�����������򣺿�ʼʱ�䡢���ʱ�䣺ƽ����ͬ�Ĳ�ֵ
                            start{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+shift_value;
                            finish{s1(Index1),s2(Index1)}=finish{s1(Index1),s2(Index1)}+shift_value;
                        end
                        %���¿���ʱ��
                        Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)}-finish{s1(LastOp),s2(LastOp)};
                    end
                    
                    %����ÿ��λ�õļӹ����ʱ��
                    Machine(mm(Index1)).MFT(j)=finish{s1(Index1),s2(Index1)};
                    if Machine(mm(Index1)).GapT(j)<0
                        Machine(mm(Index1)).GapT(j)
                    end
                end
                %���»����깤ʱ��
                mt{mm(Index1)}=Machine(mm(Index1)).MFT(ON);
            else%��290�����index1==0˵��û�п�λ��Ҫ����ʵʵ��ȥ�ں���
                %���ǻ����ϵĵ�һ���ӹ����򣬼�׼��ʱ�䣬��һ��������s1(Machine(mm(i)).Op(end))
                start{s1(i),s2(i)}=max(Machine(mm(i)).MFT(ON),finish{s1(i),s2(i)-1});
                mt{mm(i)}=start{s1(i),s2(i)}+time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(Machine(mm(i)).Op(end)),s1(i));
                finish{s1(i),s2(i)}=mt{mm(i)};
                Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
                gap=start{s1(i),s2(i)}-Machine(mm(i)).MFT(ON);
                Machine(mm(i)).GapT=[Machine(mm(i)).GapT,gap];
                Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
            end
        else%��300��ͬһ������/��ǰ������һ���ӹ����򣬲�����׼��ʱ��
            mt{mm(i)}=finish{s1(i),s2(i)-1}+time{f_index,s1(i),s2(i),mm(i)};
            start{s1(i),s2(i)}=finish{s1(i),s2(i)-1};
            finish{s1(i),s2(i)}=mt{mm(i)};
            Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
            Machine(mm(i)).GapT=[Machine(mm(i)).GapT,start{s1(i),s2(i)}];%��һ���ӹ�����ǰ�Ŀ���
            Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
        end
        
    end
end
newp=s1;
end

function [newp]=AS2FAS(p_chrom,m_chrom,FJ,f_index)%ȫ��������
global N H TM time s_time;

e=0;
JOBN=length(FJ);
SH=length(p_chrom);
finish={};%�����깤ʱ��
start={};%����ʼʱ��
for i=1:JOBN%��ʼ�����ʱ�����
    JOBI=FJ(i);
    for j=1:H(JOBI)
        finish{JOBI,j}=e;
        start{JOBI,j}=e;
    end
end

mt=cell(1,TM);
for i=1:TM%��ʼ������������ʱ������
    Machine(i).Op=[];
    Machine(i).GapT=[];
    Machine(i).MFT=[];
    mt{i}=e;
end
s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
    s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
end

for i=1:SH
    t1=s1(i);%��¼����ǰ���Ǹ�����
    t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
    mm(i)=m_chrom(1,sum(H(1,1:t1-1))+t2);%��ȡ�ù���ôμӹ��Ļ���ѡ����Ϊ����������б�ʾ�ù����ڼ��μӹ���ѡ�Ļ�������һ�α�ʾһ������
end
%��ʼ����
for i=SH:-1:1 %�������Ⱦɫ��
    if(s2(i)==H(s1(i)))
        ON=length(Machine(mm(i)).Op);%�û�����Ŀǰ�Ĺ�����
        if ON>0%355:ĳ�����ĵ�һ�����򣬵���ǰ����֮ǰ�мӹ�
            %��Ҫ�޸�t=time{f_index,s1(i),s2(i),mm(i)};�����Ǽӹ�ʱ��+׼��ʱ������бȽ�
            %t:1*ON
            t=zeros(1,ON);
            t(1)=time{f_index,s1(i),s2(i),mm(i)};%��һ������
            if ON>1
                for ON_index=2:ON
                    %��ǰ�ӹ��Ĺ���s1(Machine(mm(i)).Op(ON_index-1))
                    t(ON_index)=time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(i),s1(Machine(mm(i)).Op(ON_index-1)));
                end
            end
            %357:
            Index1=0;
            for j=1:ON%��ǰ�����ҿռ�
                if Machine(mm(i)).GapT(j)-t(j)>0
                    Index1=j; %��Index����ǰ���������Ĺ���
                    break;
                end
            end
            if Index1~=0%�ɲ��룬���º������еļӹ��Ĺ�����Ϣ
                Index1=Machine(mm(i)).Op(Index1);% ��s1�����а�s1(i)���뵽s1(Index1)����,i<Index1
                
                %����s1
                tmp=s1(i);
                for j=i:Index1-1
                    s1(j)=s1(j+1);
                end
                s1(Index1)=tmp;
                
                %����s2
                tmp=s2(i);
                for j=i:Index1-1
                    s2(j)=s2(j+1);
                end
                s2(Index1)=tmp;
                
                %���¼ӹ�����mm��ͬ�ϣ�
                tmp=mm(i);
                for j=i:Index1-1
                    mm(j)=mm(j+1);
                end
                mm(Index1)=tmp;
                
                %���µ�ǰ�����ļӹ�����
                for j=1:ON
                    if Machine(mm(Index1)).Op(j)<=Index1
                        Machine(mm(Index1)).Op(j)=Machine(mm(Index1)).Op(j)-1;%����ǰ����빤������������Ҫ���������һ��
                    end
                end
                
                %�����������л����ļӹ�����
                for k=1:TM
                    if k~=mm(Index1)
                        ON2=length(Machine(k).Op);
                        for h=1:ON2
                            if Machine(k).Op(h)<Index1&&Machine(k).Op(h)>i
                                Machine(k).Op(h)=Machine(k).Op(h)-1;
                            end
                        end
                    end
                end
                
                %���µ�ǰ�����ļӹ�������������Ĺ���
                Machine(mm(Index1)).Op=[Machine(mm(Index1)).Op,Index1];
                tmp2=Machine(mm(Index1)).Op;
                Machine(mm(Index1)).Op=sort(tmp2,'descend');
                
                %���´����빤��Ŀ�ʼ�ӹ�ʱ��
                IIndex=find(Machine(mm(Index1)).Op==Index1);
                if IIndex==1%��Ϊ��ǰ�����ĵ�һ���ӹ�����ʼʱ��Ϊ0
                    start{s1(Index1),s2(Index1)}=0;
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t(1);
                else%���¿�ʼʱ�䣨��ʼ׼������max{�ù���ǰһ������Ľ���ʱ�䣬��ǰ����ǰһ������Ľ���ʱ��}
                    LastOp=Machine(mm(Index1)).Op(IIndex-1);%�����λ��ǰһ���ӹ�����
                    start{s1(Index1),s2(Index1)}=max(0,finish{s1(LastOp),s2(LastOp)});
                    %���´����빤��Ľ����ӹ�ʱ�䣻����׼��ʱ��
                    %��һ��������s1(LastOp)
                    t_newset=time{f_index,s1(Index1),s2(Index1),mm(Index1)}+s_time(s1(Index1),s1(LastOp));%׼��ʱ��+�ӹ�ʱ��
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t_newset;
                end
                %finish{s1(Index1),s2(Index1)}=t+start{s1(Index1),s2(Index1)};
                
                %��¼���빤��
                insert_index=Index1;
                
                ON=ON+1;
                shift_value=0;
                for j=1:ON                    
                    Index1=Machine(mm(Index1)).Op(j);
                    if j==1
                        Machine(mm(Index1)).GapT(j)=0;
                    else%���ǵ�һ���ӹ�����ʱ
                        %�����ڲ��빤���ĵ�һ�����򣺿�ʼʱ�䲻�䣬׼��ʱ��ı�-->����ʱ��ı䣬ƽ��finish�Ĳ�ֵ
                        %�����ڲ��빤������������,���¹���Ŀ�ʼʱ�䡢���ʱ�䣺ƽ����ͬ�Ĳ�ֵ
                        %��������׼��ʱ�䣬���뵽����λ�ã�����Ĺ���ʼʱ�䲻��
                        LastOp=Machine(mm(Index1)).Op(j-1);%��һ������
                        if LastOp==insert_index%�����ڲ��빤���ĵ�һ������
                            new_finish=start{s1(Index1),s2(Index1)}+time{f_index,s1(Index1),s2(Index1),mm(Index1)}+s_time(s1(Index1),s1(LastOp));
                            shift_value=new_finish-finish{s1(Index1),s2(Index1)};
                            finish{s1(Index1),s2(Index1)}=new_finish;
                            %����ƽ�ƵĹ������й���s1(Index1)�Ĺ����������������ϣ�ҲҪ����start��finish��Machine.GapT��MFT
                            flag=zeros(1,TM);
                            flag(mm(Index1))=1;
                            for kk=Index1-1:-1:i
                                if s1(kk)==s1(Index1)
                                    follow_seq=kk;%O93�ļӹ�˳��=13����MAIN�������rng(0);
                                    if mm(follow_seq)~=mm(Index1)&&flag(mm(follow_seq))==0%���O93�ļӹ�������O94��ͬ���һ���û�������������O93�������ڻ�����GapT��MFT
                                        index_in_M=find(Machine(mm(follow_seq)).Op==follow_seq);
                                        for kkk=index_in_M:length(Machine(mm(follow_seq)).Op)
                                            seq=Machine(mm(follow_seq)).Op(kkk);
                                            start{s1(seq),s2(seq)}=start{s1(seq),s2(seq)}+shift_value;
                                            finish{s1(seq),s2(seq)}=finish{s1(seq),s2(seq)}+shift_value;
                                            Machine(mm(seq)).GapT(index_in_M)=Machine(mm(seq)).GapT(index_in_M)+shift_value;
                                            Machine(mm(seq)).MFT(index_in_M)=Machine(mm(seq)).MFT(index_in_M)+shift_value;
                                        end
                                        mt{mm(follow_seq)}=mt{mm(follow_seq)}+shift_value;
                                        flag(mm(follow_seq))=1;
                                    end
                                end
                            end                            
                        end
                        if LastOp<insert_index%�����ڲ��빤�����������򣺿�ʼʱ�䡢���ʱ�䣺ƽ����ͬ�Ĳ�ֵ
                            start{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+shift_value;
                            finish{s1(Index1),s2(Index1)}=finish{s1(Index1),s2(Index1)}+shift_value;
                        end
                        
                        %���¿���ʱ��
                        Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)}-finish{s1(LastOp),s2(LastOp)};
                        if Machine(mm(Index1)).GapT(j)<0
                            Machine(mm(Index1)).GapT(j)
                        end
                    end
                    %����ÿ��λ�õļӹ����ʱ��
                    Machine(mm(Index1)).MFT(j)=finish{s1(Index1),s2(Index1)};
                end
                
                %���»����깤ʱ��
                mt{mm(Index1)}=Machine(mm(Index1)).MFT(ON);
                
            else %424:���index1==0˵��û�п�λ��Ҫ����ʵʵ��ȥ�ں���
                %��һ���ӹ�������s1(Machine(mm(i)).Op(end))
                start{s1(i),s2(i)}=Machine(mm(i)).MFT(ON);%Ҳ�ǿ�ʼ׼����ʱ��
                %�����깤ʱ��+׼��ʱ��
                mt{mm(i)}=start{s1(i),s2(i)}+time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(i),s1(Machine(mm(i)).Op(end)));
                finish{s1(i),s2(i)}=mt{mm(i)};
                Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
                Machine(mm(i)).GapT=[Machine(mm(i)).GapT,0];%�޿���
                Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
            end
            
        else%433:ĳ�����ĵ�һ�������ҵ�ǰ������һ�μӹ������¹����Ŀ�ʼ������ʱ�䣬�����ӹ��Ĺ������ʱ��
            mt{mm(i)}=time{f_index,s1(i),s2(i),mm(i)};
            start{s1(i),s2(i)}=0;
            finish{s1(i),s2(i)}=mt{mm(i)};
            Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
            Machine(mm(i)).GapT=[Machine(mm(i)).GapT,0];%����ʱ��
            Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];%���ʱ��
        end
        
    else
        %443:��ʼ����ͬ�Ļ�����Ѱ�Һ��ʵĿ�λ���룬�������ļӹ�ʱ��С�ڿ���ʱ������Բ���
        ON=length(Machine(mm(i)).Op);%�û�����Ŀǰ�Ĺ�����
        if ON>0
            %��Ҫ�޸�t=time{f_index,s1(i),s2(i),mm(i)};�����Ǽӹ�ʱ��+׼��ʱ������бȽ�
            %t:1*ON
            t=zeros(1,ON);
            t(1)=time{f_index,s1(i),s2(i),mm(i)};%��һ�����У���Ϊ��һ���ӹ�������׼��ʱ��
            if ON>1
                for ON_index=2:ON
                    %��ǰ�ӹ��Ĺ���Machine(mm(i)).Op(ON_index-1)
                    t(ON_index)=time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(Machine(mm(i)).Op(ON_index-1)),s1(i));
                end
            end            
            
            Index1=0;
            for j=1:ON
                if Machine(mm(i)).GapT(j)>t(j)
                    if ON==1 || j==1
                        tmp=finish{s1(i),s2(i)+1}-0;
                    else
                        tmp=finish{s1(i),s2(i)+1}-Machine(mm(i)).MFT(j-1);
                    end
                    if Machine(mm(i)).GapT(j)-t(j)-tmp>0
                        Index1=j; %��Index����ǰ���������Ĺ���
                        break;
                    end
                end
            end
            if Index1~=0
                %462:������ĳ�����ĵ�һ�����򣬿ɲ���
                Index1=Machine(mm(i)).Op(Index1);% ��s1�����а�s1(i)���뵽s1(Index1)ǰ��,Index1<i
                
                %����s1
                tmp=s1(i);
                for j=i:Index1-1
                    s1(j)=s1(j+1);
                end
                s1(Index1)=tmp;
                
                %����s2
                tmp=s2(i);
                for j=i:Index1-1
                    s2(j)=s2(j+1);
                end
                s2(Index1)=tmp;
                
                %���¼ӹ�����mm��ͬ�ϣ�
                tmp=mm(i);
                for j=i:Index1-1
                    mm(j)=mm(j+1);
                end
                mm(Index1)=tmp;
                
                %���µ�ǰ�����ļӹ�����
                for j=1:ON
                    if Machine(mm(Index1)).Op(j)<=Index1
                        Machine(mm(Index1)).Op(j)=Machine(mm(Index1)).Op(j)-1;%����ǰ����빤������������Ҫ���������һ��
                    end
                end
                
                %�����������л����ļӹ�����
                for k=1:TM
                    if k~=mm(Index1)
                        ON2=length(Machine(k).Op);
                        for h=1:ON2
                            if Machine(k).Op(h)<Index1&&Machine(k).Op(h)>i
                                Machine(k).Op(h)=Machine(k).Op(h)-1;
                            end
                        end
                    end
                end
                
                %���µ�ǰ�����ļӹ�������������Ĺ���
                Machine(mm(Index1)).Op=[Machine(mm(Index1)).Op,Index1];
                tmp2=Machine(mm(Index1)).Op;
                Machine(mm(Index1)).Op=sort(tmp2,'descend');
                
                %���´����빤��Ŀ�ʼ�ӹ�ʱ��
                IIndex=find(Machine(mm(Index1)).Op==Index1);
                if IIndex==1%��Ϊ��ǰ�����ĵ�һ���ӹ�����ʼʱ��Ϊ0
                    start{s1(Index1),s2(Index1)}=max(0,finish{s1(Index1),s2(Index1)+1});
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t(1);
                else%���¿�ʼʱ�䣨��ʼ׼������max{�ù���ǰһ������Ľ���ʱ�䣬��ǰ����ǰһ������Ľ���ʱ��}
                    LastOp=Machine(mm(Index1)).Op(IIndex-1);
                    start{s1(Index1),s2(Index1)}=max(finish{s1(Index1),s2(Index1)+1},finish{s1(LastOp),s2(LastOp)});
                    %���´����빤��Ľ����ӹ�ʱ�䣻����׼��ʱ��
                    %��һ��������s1(LastOp)
                    t_newset=time{f_index,s1(Index1),s2(Index1),mm(Index1)}+s_time(s1(Index1),s1(LastOp));%׼��ʱ��+�ӹ�ʱ��
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t_newset;
                end
                %finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t;
                
                %��¼���빤��
                insert_index=Index1;
                
                %���µ�ǰ�����ļӹ���Ϣ
                ON=ON+1;
                shift_value=0;
                for j=1:ON                    
                    Index1=Machine(mm(Index1)).Op(j);
                    if j==1
                        Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)};
                    else%���ǵ�һ���ӹ�����ʱ
                        %�����ڲ��빤���ĵ�һ�����򣺿�ʼʱ�䲻�䣬׼��ʱ��ı�-->����ʱ��ı䣬ƽ��finish�Ĳ�ֵ
                        %�����ڲ��빤������������,���¹���Ŀ�ʼʱ�䡢���ʱ�䣺ƽ����ͬ�Ĳ�ֵ
                        %��������׼��ʱ�䣬���뵽����λ�ã�����Ĺ���ʼʱ�䲻��
                        LastOp=Machine(mm(Index1)).Op(j-1);%��һ������
                        if LastOp==insert_index%�����ڲ��빤���ĵ�һ������
                            new_finish=start{s1(Index1),s2(Index1)}+time{f_index,s1(Index1),s2(Index1),mm(Index1)}+s_time(s1(Index1),s1(LastOp));
                            shift_value=new_finish-finish{s1(Index1),s2(Index1)};
                            finish{s1(Index1),s2(Index1)}=new_finish;
                        end
                        if LastOp<insert_index%�����ڲ��빤�����������򣺿�ʼʱ�䡢���ʱ�䣺ƽ����ͬ�Ĳ�ֵ
                            start{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+shift_value;
                            finish{s1(Index1),s2(Index1)}=finish{s1(Index1),s2(Index1)}+shift_value;
                        end
                        %���¿���ʱ��
                        Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)}-finish{s1(LastOp),s2(LastOp)};
                    end
                    
                    %����ÿ��λ�õļӹ����ʱ��
                    Machine(mm(Index1)).MFT(j)=finish{s1(Index1),s2(Index1)};
                    if Machine(mm(Index1)).GapT(j)<0
                        Machine(mm(Index1)).GapT(j)
                    end
                end
                %���»����깤ʱ��
                mt{mm(Index1)}=Machine(mm(Index1)).MFT(ON);
            else%���index1==0˵��û�п�λ��Ҫ����ʵʵ��ȥ�ں���
                %���ǻ����ϵĵ�һ���ӹ����򣬼�׼��ʱ�䣬��һ��������s1(Machine(mm(i)).Op(end))
                start{s1(i),s2(i)}=max(Machine(mm(i)).MFT(ON),finish{s1(i),s2(i)+1});
                mt{mm(i)}=start{s1(i),s2(i)}+time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(i),s1(Machine(mm(i)).Op(end)));
                finish{s1(i),s2(i)}=mt{mm(i)};
                Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
                gap=start{s1(i),s2(i)}-Machine(mm(i)).MFT(ON);
                Machine(mm(i)).GapT=[Machine(mm(i)).GapT,gap];
                Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
            end
        else%531��ͬһ������/��ǰ������һ���ӹ����򣬲�����׼��ʱ��
            mt{mm(i)}=finish{s1(i),s2(i)+1}+time{f_index,s1(i),s2(i),mm(i)};
            start{s1(i),s2(i)}=finish{s1(i),s2(i)+1};
            finish{s1(i),s2(i)}=mt{mm(i)};
            Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
            Machine(mm(i)).GapT=[Machine(mm(i)).GapT,start{s1(i),s2(i)}];
            Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
        end
        
    end
end

newp=s1;
end