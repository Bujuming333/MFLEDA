function [fit1,fit2,fit3]=fitDFJSP(p_chrom,m_chrom,f_chrom)
%[fit1,fit2,fit3]������깤ʱ�䣬���ܺģ��ؼ�����
global N SH F;
%fit1:Cmax
%fit2�����ܺ�
%fit3���ؼ�����[0��F-1]
s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
    s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
end

P=cell(1,F);%P=1*3 cell��=[] [] [],ÿ�������ļӹ�����
for i=1:SH%ÿ�������ӹ��Ĺ����ţ���initialʱ�Ѿ��̶��ӹ�����������6,7,7��P=[1*30][1*35][1*35]
    t1=s1(i);%��¼����ǰ���Ǹ�����
    t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
    P{f_chrom(t1)+1}=[P{f_chrom(t1)+1} p_chrom(i)];
end
FJ=cell(1,F);%FJÿ�������ļӹ�����
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end
Fit=zeros(F,2);%Fit��F*2��ÿ������������Ŀ��ֵ
for i=1:F
    [Fit(i,1),Fit(i,2)]=fitFJSP(P{i},m_chrom,FJ{i},i);
end

%��ע˭�ǹؼ�����
[~,fit3]=max(Fit(:,1));
fit3=fit3-1;

fit1=max(Fit(:,1));
fit2=sum(Fit(:,2));
end

function [fit1,fit2]=fitFJSP(pro_m,mac_m,FJ,F_index)%���ô�ͳ��˫����뷽ʽ ������깤ʱ��makespan ��ʾ���һ����ɹ����ʱ��
global N H  TM time s_time;

Ep=4;%�ӹ����� 500kwÿs 0.5Mw
Ei=1;%�ȴ����� 40kwÿs
Es=2;%׼��ʱ�书��
%     Ec=0.5;%���������ȹ��õĹ��� 60kwÿs

JOBN=length(FJ);
SH=length(pro_m);

finish={};%�����깤ʱ��
for i=1:JOBN%��ʼ�����ʱ�����
    JOBI=FJ(i);
    for j=1:H(JOBI)
        finish{i,j}=0;
    end
end

mt=cell(1,TM);
for i=1:TM%��ʼ������������ʱ������
    mt{i}=0;
end

s1=pro_m;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
    s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
end

for i=1:SH
    t1=s1(i);%��¼����ǰ���Ǹ�����
    t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
    mm(i)=mac_m(1,sum(H(1,1:t1-1))+t2);%��ȡ�ù���ôμӹ��Ļ���ѡ����Ϊ����������б�ʾ�ù����ڼ��μӹ���ѡ�Ļ�������һ�α�ʾһ������
end
total=0;%�ӹ�ʱ��
I_total=0;%�ܵĹ���ȴ�ʱ��
s_total=0;%׼��ʱ��
for i=1:SH
    if(s2(i)==1)
        %�ҵ�ǰ������ǰһ���ӹ�����
        pre_i=i;
        if i>1
            for find_i=i-1:-1:1
                if mm(find_i)==mm(i)
                    pre_i=find_i;
                    break;
                end
            end
        end
        
        mt{mm(i)}=mt{mm(i)}+time{F_index,s1(i),s2(i),mm(i)}+s_time(s1(pre_i),s1(i));%�ۼƼ���ÿһ̨�����ļӹ�ʱ��
        finish{s1(i),s2(i)}= mt{mm(i)};%�Ĺ�����ǰ��������ʱ���ǵ�ǰ�������ۼ�ʱ��
        total=total+time{F_index,s1(i),s2(i),mm(i)};%�ܵĻ�������
        s_total=s_total+s_time(s1(pre_i),s1(i));
    else
        if(mt{mm(i)}<finish{s1(i),s2(i)-1})%�������������ʱ��С�ڸù�����һ�����ʱ��
            
            %�ҵ�ǰ������ǰһ���ӹ�����
            pre_i=i;
            if i>1
                for find_i=i-1:-1:1
                    if mm(find_i)==mm(i)
                        pre_i=find_i;
                        break;
                    end
                end
            end
            Idletime=finish{s1(i),s2(i)-1}-mt{mm(i)};
            I_total=I_total+Idletime;%���ʱ��������ֿ���ʱ�䣬����Ҫ���Ĺ��ʡ�
            
            mt{mm(i)}= finish{s1(i),s2(i)-1}+time{F_index,s1(i),s2(i),mm(i)}+s_time(s1(pre_i),s1(i));%һ����ȥ��һ�����ʱ��͸û��������ӹ�ʱ��������
            finish{s1(i),s2(i)}= mt{mm(i)};
            total=total+time{F_index,s1(i),s2(i),mm(i)};
            s_total=s_total+s_time(s1(pre_i),s1(i));
        else%����������ʱ����ڵ��ڸù�����һ����������ʱ��
            
            %�ҵ�ǰ������ǰһ���ӹ�����
            pre_i=i;
            if i>1
                for find_i=i-1:-1:1
                    if mm(find_i)==mm(i)
                        pre_i=find_i;
                        break;
                    end
                end
            end
            
            mt{mm(i)}= mt{mm(i)}+time{F_index,s1(i),s2(i),mm(i)}+s_time(s1(pre_i),s1(i));
            finish{s1(i),s2(i)}= mt{mm(i)};
            total=total+time{F_index,s1(i),s2(i),mm(i)};
            s_total=s_total+s_time(s1(pre_i),s1(i));
        end
    end
end

fit1=mt{1};
for i=2:TM
    if(mt{i}>fit1)%�������i������깤ʱ�����fit�����
        fit1=mt{i};
    end
end
fit2=total*Ep+I_total*Ei+s_total*Es; %�ܺ�

end