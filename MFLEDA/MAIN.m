%数据集：DHFJSP; SDST
%distributed heterogeneous flexible job shop scheduling problem with sequence-dependent setup times
%the same operation processed in different factory has the same candidate machine selection but different processing times
%the whole operations of a job must be processed on the same factory
%测试集DHFJSP来源：Li Rui-2023-Co-evolution with Deep Reinforcement Learning For Energy-Aware Distributed Heterogeneous Flexible Job Shop Scheduling
%三种不同的方式，UMDA;PBIL;GMEDA更新工厂分配，机器选择，工序加工顺序
%局部多目标地形特征自适应选择算子
clc;
clear;
DataPath='DATASET\';
job_type=[10,20,30,40,50,100,150,200];
f_type={[2],[2,3],[2,3],[2,3,4],[3,4,5],[4,5,6,7],[5,6,7],[6,7]};
nameJ='J';nameF='F'; txt='.txt';
paint=1;%%画图flag
k=1;
for i=1:size(job_type,2)
    for j=f_type{i}
        RealPath{k}=[DataPath,num2str(job_type(i)),nameJ,num2str(j),nameF,txt];
        ResultPath{k}=[num2str(job_type(i)),nameJ,num2str(j),nameF];
        k=k+1;
    end
end

global  N H SH NM ps M TM time AP AM AF AN F ProF s_time;

RP=[[0,0],[5,0],[10,0],[0,5],[0,10]];%GD指标需要的参考点

for File=1:20%20个实例
    %%% NM：N*5，1~5，可加工工序O_ij的机器数，最大为5
    %%% M：N*TM*5，1~5，可加工工序O_ij的机器号
    [N,F,TM,H,SH,NM,M,time,ProF]=DataReadDHFJSP(RealPath{File});
    MaxNFEs=200*SH;
    ps=100;%{100,150,200}
    alpha_max=0.5;%{0.2,0.3,0.4,0.5}
    length_l=5;%{20,50,100,150,200}随机游走长度
    
    respath='result\';
    tmp6='\';
    respath=[respath,ResultPath{File},tmp6];
    %清除无关变量
    clear tmp6
    a=['mkdir ' respath];%创立写结果的文件夹
    system(a);%利用Windows命令行命令执行dos命令
    fprintf('%s %s\r\n','Calculating ',RealPath{File});
    
    %读取准备时间
    FileName1=['SDST\SDST_',strcat(int2str(N)),'.txt'];
    s_time=textread(FileName1);
    
    for round_id=1:20
        %初始化
        [p_chrom,m_chrom,f_chrom]=initial();
        %评价种群
        NFEs=0;
        fitness=zeros(ps,3);
        for i=1:ps
            [fitness(i,1),fitness(i,2),fitness(i,3)]=fitDFJSP(p_chrom(i,:),m_chrom(i,:),f_chrom(i,:));
            NFEs=NFEs+1;
        end
        AP=[];AM=[];AF=[];AN=[];
        A=cell(1,N);%Oi,j分配到第k个可选机器
        for i=1:N
            A{i}=zeros(H(i),5);
            for j=1:H(i)
                for k=1:NM{i,j}
                    A{i}(j,k)=1/NM{i,j};
                end
            end
        end
        iter_index=1;
        while NFEs<MaxNFEs %最大评价次数
            fprintf('%s %s %d %s %d\r\n',RealPath{File},'round',round_id,'iter',iter_index);
            iter_index=iter_index+1;
            
            %选择——随机选择两个比较支配关系，所以有重复
            [Par_p,Par_m,Par_f,Par_Fit]=TSelection(p_chrom,m_chrom,f_chrom,fitness);
            CP=[];CM=[];CN=[];CF=[];
            
            %建模、采样
            %%%%%%GMEDA
            [new_p1]=GMEDA_p(Par_p);
            %%%%%%PBIL
            alpha=alpha_max-NFEs/MaxNFEs*(alpha_max-0.01);
            [new_m1,A]=PBIL_m(Par_m,A,alpha);
            %%%%%%UMDA
            [new_f1]=UMDA_f(Par_f);%与initial不同，每个解中工厂分配到的工件数不再均匀
            
            %评价
            CF=zeros(ps,3);
            for j=1:ps
                [CF(j,1),CF(j,2),CF(j,3)]=fitDFJSP(new_p1(j,:),new_m1(j,:),new_f1(j,:));
                NFEs=NFEs+1;
            end
            CP=new_p1;CM=new_m1;CN=new_f1;
            
            %环境选择
            QP=[];QM=[];QN=[];QF=[];
            QP=[Par_p;CP];QM=[Par_m;CM];QN=[Par_f;CN];QF=[Par_Fit;CF];%文章——合并Q_t（交配池）和C_t；ps+2*ps,共3*ps个个体；修改后选择的ps个个体+新生成的ps个个体，共2*ps个个体
            [QP,QM,QN,QF]=DeleteRpeatQF(QP,QM,QN,QF);
            [TopPSRank]=FastNDS(QF);
            
            %更新新种群
            p_chrom=QP(TopPSRank,:);m_chrom=QM(TopPSRank,:);f_chrom=QN(TopPSRank,:);fitness=QF(TopPSRank,:);
            for j=1:ps%更新AP,AM,AN,AF；全主动调度的能量存储策略，并将更新了父代的子代解加入存档Ω
                [p_chrom(j,:),m_chrom(j,:),f_chrom(j,:),fitness(j,:)]=EnergySaveDFJSP(p_chrom(j,:),m_chrom(j,:),f_chrom(j,:),fitness(j,:));
                NFEs=NFEs+1;
            end
            
            %Update Elite Archive
            [PF,~]=pareto1(fitness);
            AP=[AP;p_chrom(PF,:)];AM=[AM;m_chrom(PF,:)];AN=[AN;f_chrom(PF,:)];AF=[AF;fitness(PF,:)];
            [PF,~]=pareto1(AF);
            AP=AP(PF,:);AM=AM(PF,:);AN=AN(PF,:);AF=AF(PF,:);
            DeleteRpeat();
            
            %局部搜索
            numst=6;
            [L,~]=size(AP);
            %Local multiobjective landscape features
            prob=RW(length_l,numst);
            NFEs=NFEs+numst*length_l;
            for j=1:L
                %根据概率prob轮盘赌选择
                chooseA=0;
                r=rand;
                while r>prob(chooseA+1)
                    chooseA=chooseA+1;
                end
                chooseA=chooseA+1;
                switch chooseA
                    case{1}
                        [P1,M1,N1]=LS1(AP(j,:),AM(j,:),AN(j,:),AF(j,:));%生成邻域解
                    case{2}
                        [P1,M1,N1]=LS2(AP(j,:),AM(j,:),AN(j,:),AF(j,:));%生成邻域解
                    case{3}
                        [P1,M1,N1]=LS3(AP(j,:),AM(j,:),AN(j,:),AF(j,:));%生成邻域解
                    case{4}
                        [P1,M1,N1]=LS4(AP(j,:),AM(j,:),AN(j,:),AF(j,:));%生成邻域解
                    case{5}
                        [P1,M1,N1]=LS5(AP(j,:),AM(j,:),AN(j,:),AF(j,:));%生成邻域解
                    case{6}
                        [P1,M1,N1]=LS6(AP(j,:),AM(j,:),AN(j,:),AF(j,:));%生成邻域解
                end
                [F1(1,1),F1(1,2),F1(1,3)]=fitDFJSP(P1,M1,N1);
                if NDS(F1,AF(j,:))==1 %新解支配旧解
                    AP(j,:)=P1;
                    AM(j,:)=M1;
                    AN(j,:)=N1;
                    AF(j,:)=F1;
                    AP=[AP;P1];
                    AM=[AM;M1];
                    AF=[AF;F1];
                    AN=[AN;N1];
                elseif NDS(F1,AF(j,:))==0 %新解互不支配旧解
                    if F1(1,1)~=AF(j,1)&&F1(1,2)~=AF(j,2)
                        AP=[AP;P1];
                        AM=[AM;M1];
                        AF=[AF;F1];
                        AN=[AN;N1];
                    end
                end
            end %JEND
        end%END WHILE
        
        %更新存档
        [PF,~]=pareto1(AF);
        AP=AP(PF,:);AM=AM(PF,:);AN=AN(PF,:);AF=AF(PF,:);
        %%添加：画甘特图
        my_drawDFJSP(AP(1,:),AM(1,:),AN(1,:));
        
        [PF,~]=pareto1(AF);
        L=length(PF);
        obj=AF(:,1:2);
        newobj=[];
        for i=1:L
            newobj(i,:)=obj(PF(i),:);
        end
        newobj=unique(newobj,'rows');%删除重复行
        tmp5=newobj';%将求得的前沿解集转置
        %scatter(newobj(:,1),newobj(:,2));
        tmp1='res';%根据字符串编写出结果存储的路径
        resPATH=[respath tmp1 num2str(round_id) txt];
        fout=fopen(resPATH,'w');
        fprintf(fout,'%5.2f %6.3f\r\n',tmp5);%由于matlab是按照列存储和读取矩阵，则写入文件的时候也是按照列填写，因此要转置存入
        fclose(fout);
    end
end