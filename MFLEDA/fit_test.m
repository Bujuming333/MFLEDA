function []=fit_test(p_chrom,m_chrom,f_chrom)
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

P=cell(1,F);
for i=1:SH
    t1=s1(i);%记录到当前是那个工件
    t2=s2(i);%记录当前工件是加工到第几次
    P{f_chrom(t1)+1}=[P{f_chrom(t1)+1} p_chrom(i)];
end
FJ=cell(1,F);
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end

for i=1:F
    fitFJSP(P{i},m_chrom,FJ{i},i);
end



end

function fitFJSP(pro_m,mac_m,FJ,F_index)%采用传统的双层解码方式 求最大完工时间makespan 表示最后一个完成工序的时间
global N H  TM time s_time;

Ep=4;%加工功率 500kw每s 0.5Mw
Ei=1;%等待功率 40kw每s
Es=2;
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
%% 判断

for i=1:SH
    if isempty(time{F_index,s1(i),s2(i),mm(i)})
        fprintf('\n%d\t%d\t%d\t%d\n%s\n',F_index,s1(i),s2(i),mm(i),'WRONG!');
    end
end
%% 



end