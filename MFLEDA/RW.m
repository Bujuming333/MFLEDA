function prob=RW(length_l,numst)
global AP AM AN AF
%length_l 随机游走长度
%numst 邻域搜索个数
%随机选择一个解
sup_avg_rws=[];
hv_avg_rws=[];
PF=cell(1,numst);
pf_sample=cell(1,numst);
for num=1:numst
    switch num
        case{1}
            fhd=@LS1;
        case{2}
            fhd=@LS2;
        case{3}
            fhd=@LS3;
        case{4}
            fhd=@LS4;
        case{5}
            fhd=@LS5;
        case{6}
            fhd=@LS6;
    end
    P_l=size(AP,1);
    Index=ceil(rand*P_l);
    
    p_sample=repmat(AP(Index,:),length_l+1,1);
    m_sample=repmat(AM(Index,:),length_l+1,1);
    f_sample=repmat(AN(Index,:),length_l+1,1);
    fit_sample=repmat(AF(Index,:),length_l+1,1);
    num_NDS=zeros(1,length_l);
    
    for j=1:length_l
        [P1,M1,F1]=feval(fhd,p_sample(j,:),m_sample(j,:),f_sample(j,:),fit_sample(j,:));%生成邻域解
        %评价
        [Fit1(1,1),Fit1(1,2),Fit1(1,3)]=fitDFJSP(P1,M1,F1);
        p_sample(j+1,:)=P1;
        m_sample(j+1,:)=M1;
        f_sample(j+1,:)=F1;
        fit_sample(j+1,:)=Fit1;
    end
    
    %#sup_avg_rws
    for j=1:length_l
        if NDS(fit_sample(j+1,:),fit_sample(1,:))~=2 %新解支配旧解；互不支配
            num_NDS(j)=1;
        end
    end
    sup_avg_rws=[sup_avg_rws sum(num_NDS)/length_l];
    
    %hv_avg_rws
    sample=fit_sample(2:end,1:2);
    PF{num}=pareto1(sample);
    pf_sample{num}=sample(PF{num},:);
end
setAll=[];
for num=1:numst
    setAll=[setAll;pf_sample{num}];
end
MIN=min(setAll);
MAX=max(setAll);
%归一化
newPF=cell(1,numst);
for num=1:numst
    for j=1:size(pf_sample{num},1)
        temp=zeros(1,2);
        temp(1)=(pf_sample{num}(j,1)-MIN(1))/(MAX(1)-MIN(1));
        temp(2)=(pf_sample{num}(j,2)-MIN(2))/(MAX(2)-MIN(2));
        newPF{num}=[newPF{num};temp];
    end
end
optimum=[1.2,1.2];
for num=1:numst
    hv=HV(newPF{num},optimum);
    hv_avg_rws=[hv_avg_rws hv];
end

if sum(sup_avg_rws)==0
    pro_sup=ones(1,numst)*1/numst;
else
    pro_sup=sup_avg_rws/sum(sup_avg_rws);
end
pro_hv=hv_avg_rws/sum(hv_avg_rws);
pro=pro_sup+pro_hv;
prob=cumsum(pro)/sum(pro);