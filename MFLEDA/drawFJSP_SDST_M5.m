function drawFJSP_SDST_M5(pro_m,mac_m,FJ,f_index)%采用传统的双层解码方式 求最大完工时间makespan 表示最后一个完成工序的时间
%5台机器
global N H  TM time F;
FileName1=['SDST\SDST_',strcat(int2str(N)),'.txt'];
s_time=textread(FileName1);
figure
e=0;
number_size=10;
O_size=12;
dot_size=10;
O_up=-2.0;
dot_up=-2.3;
O_down=1.4;
dot_down=1.7;
mshift=-3;%-15;
%title(['Factory',strcat(int2str(f_index)),'/',strcat(int2str(F)),' of ',strcat(int2str(length(FJ))),'/',strcat(int2str(N)),' jobs']);
%title(['Factory',strcat(int2str(f_index))]);
xlen=120;
axis([0,xlen,-2,14]);
set(gca,'FontSize',12,'Fontname','Times New Roman');
hold on;
% plot([0 xlen],[2 2]);%画横线 x=0-50 y=2-2
% hold on;
% plot([0 xlen],[6 6]);
% hold on;
% plot([0 xlen],[10 10]);
% hold on;
% plot([0 xlen],[14 14]);
% hold on;
% plot([0 xlen],[18 18]);

%       hold on;
%      plot([0 70],[26 26]);
%       hold on;
%      plot([0 70],[30 30]);
%       hold on;
%      plot([0 70],[34 34]);
%       hold on;
%      plot([0 70],[38 38]);
hold on;
text(mshift,0,'M_{5}','FontSize',12,'FontName','Times New Roman');
text(mshift,3,'M_{4}','FontSize',12,'FontName','Times New Roman');
text(mshift,6,'M_{3}','FontSize',12,'FontName','Times New Roman');
text(mshift,9,'M_{2}','FontSize',12,'FontName','Times New Roman');
text(mshift,12,'M_{1}','FontSize',12,'FontName','Times New Roman');
text(-6,5,strcat('F_{',int2str(f_index),'}'),'FontSize',12,'FontName','Times New Roman');

%       text(-3,26,'M_{4}');
%      text(-3,30,'M_{3}');
%      text(-3,34,'M_{2}');
%       text(-3,38,'M_{1}')
c={[1 0.93725  0.83529] [1 0.8549   0.72549] [0.6902 0.87843 0.90196] [0.25098 0.87843 0.81569] [0 0.74902 1] [1 0.92157 0.80392] [0.73725 0.82353 0.93333] [0.93333 0.9098 0.66667] [0.52941 0.80784 0.92157] [0.93333 0.81176 0.63137]...
    [0.39216 0.58431 0.92941] [0.80392 0.70196 0.5451] [0.4 0.6 0.8] [0.11765 0.56471 1] [0.4 0.80392 0.66667] [0.96078 0.87059 0.70196] [0.28235 0.81961 0.8] [0.27451 0.5098 0.70588] [0.52941 0.80784 0.98039] [0.2549 0.41176 0.88235]};
SH=length(pro_m);
JOBN=length(FJ);

finish={};
start={};
for i=1:JOBN%初始化完成时间矩阵
    for j=1:H(i)
        finish{i,j}=e;
        start{i,j}=e;
    end
end

mt=cell(1,TM);
for i=1:TM%初始化机器最大完成时间数组
    mt{i}=e;
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
for i=1:SH
    if mm(i)==6
        mm(i)=randperm(5,1);
        if isempty(time{f_index,s1(i),s2(i),mm(i)})
            time{f_index,s1(i),s2(i),mm(i)}=randperm(6,1);
        end
    end
end

for i=1:SH
    if(s2(i)==1)
        mx=mm(i);
        mx=16-mx*3;
        
        pre_i=i;
        if i>1
            for find_i=i-1:-1:1
                if mm(find_i)==mm(i)
                    pre_i=find_i;
                    break;
                end
            end
        end
        start_setup=mt{mm(i)};
        start{s1(i),s2(i)}= mt{mm(i)}+s_time(s1(pre_i),s1(i));
        %start{s1(i),s2(i)}= mt{mm(i)};
        
        mt{mm(i)}=start{s1(i),s2(i)}+time{f_index,s1(i),s2(i),mm(i)};
        %mt{mm(i)}=mt{mm(i)}+time{f_index,s1(i),s2(i),mm(i)};%累计计算每一台机器的加工时间
        
        
        finish{s1(i),s2(i)}= mt{mm(i)};%改工件当前工序的完成时间是当前机器的累计时间
        
        %画准备时间
        if s_time(s1(pre_i),s1(i))>0
            rectangle('position',[start_setup,mx-1.5,s_time(s1(pre_i),s1(i)),1.5],'Facecolor',[0.8,0.8,0.8]);
        end
        
        s=start{s1(i),s2(i)};f=finish{s1(i),s2(i)};p=time{f_index,s1(i),s2(i),mm(i)};
        rectangle('position',[s,mx-1.5,p,1.5],'Facecolor',c{s1(i)});
        text(start{s1(i),s2(i)}+p/2-0.8,mx-0.8,1,strcat('O_{',int2str(s1(i)),strcat(',',int2str(s2(i))),'}'),'FontSize',O_size,'Fontname','Times New Roman');

    else
        if(mt{mm(i)}<finish{s1(i),s2(i)-1})%如果机器最大完成时间小于该工序上一步完成时间
            mx=mm(i);
            mx=16-mx*3;
            
            %画空闲时间
            Idletime=finish{s1(i),s2(i)-1}-mt{mm(i)};
            start_idle=mt{mm(i)};
            if mt{mm(i)}==0
                Idletime=0;
            end
            rectangle('position',[start_idle,mx-1.5,Idletime,1.5],'Facecolor',[1,1,1]);
                        
            start{s1(i),s2(i)}= finish{s1(i),s2(i)-1};
            start_setup=start{s1(i),s2(i)};%开始准备时间
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
            start{s1(i),s2(i)}=start{s1(i),s2(i)}+s_time(s1(pre_i),s1(i));
            mt{mm(i)}= start{s1(i),s2(i)}+time{f_index,s1(i),s2(i),mm(i)};
            %mt{mm(i)}= finish{s1(i),s2(i)-1}+time{f_index,s1(i),s2(i),mm(i)};%一定是去上一步完成时间和该机器结束加工时间的最大者
            finish{s1(i),s2(i)}= mt{mm(i)};
            
            %画准备时间
            if s_time(s1(pre_i),s1(i))>0
                rectangle('position',[start_setup,mx-1.5,s_time(s1(pre_i),s1(i)),1.5],'Facecolor',[0.8,0.8,0.8]);
            end
            
            s=start{s1(i),s2(i)};f=finish{s1(i),s2(i)};p=time{f_index,s1(i),s2(i),mm(i)};
            rectangle('position',[s,mx-1.5,p,1.5],'Facecolor',c{s1(i)});
            text(start{s1(i),s2(i)}+p/2-0.8,mx-0.8,1,strcat('O_{',int2str(s1(i)),strcat(',',int2str(s2(i))),'}'),'FontSize',O_size,'Fontname','Times New Roman');

        else%如果机器完成时间大于等于该工件上一个工序的完成时间
            mx=mm(i);
            mx=16-mx*3;
            start{s1(i),s2(i)}= mt{mm(i)};
            start_setup=start{s1(i),s2(i)};
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
            start{s1(i),s2(i)}=start{s1(i),s2(i)}+s_time(s1(pre_i),s1(i));
            
            mt{mm(i)}= start{s1(i),s2(i)}+time{f_index,s1(i),s2(i),mm(i)};
            %mt{mm(i)}= mt{mm(i)}+time{f_index,s1(i),s2(i),mm(i)};
            finish{s1(i),s2(i)}= mt{mm(i)};
            
            %画准备时间
            if s_time(s1(pre_i),s1(i))>0
                rectangle('position',[start_setup,mx-1.5,s_time(s1(pre_i),s1(i)),1.5],'Facecolor',[0.8,0.8,0.8]);
            end
            
            s=start{s1(i),s2(i)};f=finish{s1(i),s2(i)};p=time{f_index,s1(i),s2(i),mm(i)};
            rectangle('position',[s,mx-1.5,p,1.5],'Facecolor',c{s1(i)});
            text(start{s1(i),s2(i)}+p/2-0.8,mx-0.8,1,strcat('O_{',int2str(s1(i)),strcat(',',int2str(s2(i))),'}'),'FontSize',O_size,'Fontname','Times New Roman');

        end
    end
end

fit1=mt{1};
for i=2:TM
    if(mt{i}>fit1)%如果机器i的最大完工时间大于fit则更新
        fit1=mt{i};
    end
end
%set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w');
% xlabel(['makespan=',strcat(int2str(fit1))]);
set(gca,'ytick',[]);
alpha(0.6)
hold off
end