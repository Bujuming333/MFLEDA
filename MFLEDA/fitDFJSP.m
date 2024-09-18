function [fit1,fit2,fit3]=fitDFJSP(p_chrom,m_chrom,f_chrom)
%[fit1,fit2,fit3]：最大完工时间，总能耗，关键工厂
global N SH F;
%fit1:Cmax
%fit2：总能耗
%fit3：关键工厂[0，F-1]
s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
    s2(i)=p(s1(i));%记录加工过程中，工件的次数
end

P=cell(1,F);%P=1*3 cell，=[] [] [],每个工厂的加工工序
for i=1:SH%每个工厂加工的工件号（在initial时已经固定加工的数量）：6,7,7，P=[1*30][1*35][1*35]
    t1=s1(i);%记录到当前是那个工件
    t2=s2(i);%记录当前工件是加工到第几次
    P{f_chrom(t1)+1}=[P{f_chrom(t1)+1} p_chrom(i)];
end
FJ=cell(1,F);%FJ每个工厂的加工工件
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end
Fit=zeros(F,2);%Fit：F*2，每个工厂的两个目标值
for i=1:F
    [Fit(i,1),Fit(i,2)]=fitFJSP(P{i},m_chrom,FJ{i},i);
end

%标注谁是关键工厂
[~,fit3]=max(Fit(:,1));
fit3=fit3-1;

fit1=max(Fit(:,1));
fit2=sum(Fit(:,2));
end

function [fit1,fit2]=fitFJSP(pro_m,mac_m,FJ,F_index)%采用传统的双层解码方式 求最大完工时间makespan 表示最后一个完成工序的时间
global N H  TM time s_time;

Ep=4;%加工功率 500kw每s 0.5Mw
Ei=1;%等待功率 40kw每s
Es=2;%准备时间功率
%     Ec=0.5;%工厂照明等公用的功率 60kw每s

JOBN=length(FJ);
SH=length(pro_m);

finish={};%工序完工时间
for i=1:JOBN%初始化完成时间矩阵
    JOBI=FJ(i);
    for j=1:H(JOBI)
        finish{i,j}=0;
    end
end

mt=cell(1,TM);
for i=1:TM%初始化机器最大完成时间数组
    mt{i}=0;
end

s1=pro_m;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
    s2(i)=p(s1(i));%记录加工过程中，工件的次数
end

for i=1:SH
    t1=s1(i);%记录到当前是那个工件
    t2=s2(i);%记录当前工件是加工到第几次
    mm(i)=mac_m(1,sum(H(1,1:t1-1))+t2);%提取该工序该次加工的机器选择，因为机器码的排列表示该工件第几次加工所选的机器，是一段表示一个工件
end
total=0;%加工时间
I_total=0;%总的工序等待时间
s_total=0;%准备时间
for i=1:SH
    if(s2(i)==1)
        %找当前机器的前一个加工工件
        pre_i=i;
        if i>1
            for find_i=i-1:-1:1
                if mm(find_i)==mm(i)
                    pre_i=find_i;
                    break;
                end
            end
        end
        
        mt{mm(i)}=mt{mm(i)}+time{F_index,s1(i),s2(i),mm(i)}+s_time(s1(pre_i),s1(i));%累计计算每一台机器的加工时间
        finish{s1(i),s2(i)}= mt{mm(i)};%改工件当前工序的完成时间是当前机器的累计时间
        total=total+time{F_index,s1(i),s2(i),mm(i)};%总的机器负载
        s_total=s_total+s_time(s1(pre_i),s1(i));
    else
        if(mt{mm(i)}<finish{s1(i),s2(i)-1})%如果机器最大完成时间小于该工序上一步完成时间
            
            %找当前机器的前一个加工工件
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
            I_total=I_total+Idletime;%这个时候机器出现空闲时间，则需要消耗功率。
            
            mt{mm(i)}= finish{s1(i),s2(i)-1}+time{F_index,s1(i),s2(i),mm(i)}+s_time(s1(pre_i),s1(i));%一定是去上一步完成时间和该机器结束加工时间的最大者
            finish{s1(i),s2(i)}= mt{mm(i)};
            total=total+time{F_index,s1(i),s2(i),mm(i)};
            s_total=s_total+s_time(s1(pre_i),s1(i));
        else%如果机器完成时间大于等于该工件上一个工序的完成时间
            
            %找当前机器的前一个加工工件
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
    if(mt{i}>fit1)%如果机器i的最大完工时间大于fit则更新
        fit1=mt{i};
    end
end
fit2=total*Ep+I_total*Ei+s_total*Es; %能耗

end