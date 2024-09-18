function [newp]=compare_SAS2AS(p_chrom,m_chrom,FJ,f_index)%主动解码
global N H TM time s_time;

e=0;
JOBN=length(FJ);
SH=length(p_chrom);
finish={};%工序完工时间
start={};%工序开始时间
for i=1:JOBN%初始化完成时间矩阵
    JOBI=FJ(i);
    for j=1:H(JOBI)
        finish{JOBI,j}=e;
        start{JOBI,j}=e;
    end
end

mt=cell(1,TM);
for i=1:TM%初始化机器最大完成时间数组
    Machine(i).Op=[];
    Machine(i).GapT=[];
    Machine(i).MFT=[];
    mt{i}=e;
end
s1=p_chrom;
s2=zeros(1,SH);
p=zeros(1,N);

for i=1:SH
    p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
    s2(i)=p(s1(i));%记录加工过程中，工件的次数
end

for i=1:SH
    t1=s1(i);%记录到当前是那个工件
    t2=s2(i);%记录当前工件是加工到第几次
    mm(i)=m_chrom(1,sum(H(1,1:t1-1))+t2);%提取该工序该次加工的机器选择，因为机器码的排列表示该工件第几次加工所选的机器，是一段表示一个工件
end
%开始解码
for i=1:SH
    if(s2(i)==1)
        ON=length(Machine(mm(i)).Op);%该机器上目前的工序数
        if ON>0%某工件的第一个工序，但当前机器之前有加工
            %需要修改t=time{f_index,s1(i),s2(i),mm(i)};：考虑加工时间+准备时间与空闲比较
            %t:1*ON
            t=zeros(1,ON);
            t(1)=time{f_index,s1(i),s2(i),mm(i)};%第一个空闲          
            if ON>1
                for ON_index=2:ON
                    %当前加工的工序s1(Machine(mm(i)).Op(ON_index-1))
                    t(ON_index)=time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(Machine(mm(i)).Op(ON_index-1)),s1(i));
                end
            end
            %127：
            Index1=0;
            for j=1:ON%从前到后找空间
                if Machine(mm(i)).GapT(j)-t(j)>0
                    Index1=j; %在Index工序前插入新来的工序
                    break;
                end
            end
            if Index1~=0%可插入，更新后面所有的加工的工序信息
                Index1=Machine(mm(i)).Op(Index1);% 在s1向量中把s1(i)插入到s1(Index1)前面,Index1<i
                
                %更新s1
                tmp=s1(i);
                for j=i:-1:Index1+1
                    s1(j)=s1(j-1);
                end
                s1(Index1)=tmp;
                
                %更新s2
                tmp=s2(i);
                for j=i:-1:Index1+1
                    s2(j)=s2(j-1);
                end
                s2(Index1)=tmp;
                
                %更新加工机器mm（同上）
                tmp=mm(i);%当前机器
                for j=i:-1:Index1+1
                    mm(j)=mm(j-1);
                end
                mm(Index1)=tmp;
                
                %更新当前机器的加工索引
                for j=1:ON
                    if Machine(mm(Index1)).Op(j)>=Index1
                        Machine(mm(Index1)).Op(j)=Machine(mm(Index1)).Op(j)+1;%由于前面插入工序，其他索引需要整体向后退一格
                    end
                end
                
                %更新其他所有机器的加工索引
                for k=1:TM
                    if k~=mm(Index1)
                        ON2=length(Machine(k).Op);%此刻机器k加工的总次数
                        for h=1:ON2
                            if Machine(k).Op(h)>Index1&&Machine(k).Op(h)<i
                                Machine(k).Op(h)=Machine(k).Op(h)+1;
                            end
                        end
                    end
                end
                
                %更新当前机器的加工索引（待插入的工序）
                Machine(mm(Index1)).Op=[Machine(mm(Index1)).Op,Index1];
                tmp2=Machine(mm(Index1)).Op;
                Machine(mm(Index1)).Op=sort(tmp2,'ascend');
                
                %更新待插入工序的开始加工时间
                IIndex=find(Machine(mm(Index1)).Op==Index1);
                if IIndex==1%若为当前机器的第一个加工，开始时间为0
                    start{s1(Index1),s2(Index1)}=0;
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t(1);
                else%更新开始时间（开始准备）：max{该工件前一个工序的结束时间，当前机器前一个工序的结束时间}
                    LastOp=Machine(mm(Index1)).Op(IIndex-1);%插入的位置前一个加工工序
                    start{s1(Index1),s2(Index1)}=max(0,finish{s1(LastOp),s2(LastOp)});
                    %更新待插入工序的结束加工时间；考虑准备时间
                    %上一个工件：s1(LastOp)
                    t_newset=time{f_index,s1(Index1),s2(Index1),mm(Index1)}+s_time(s1(LastOp),s1(Index1));%准备时间+加工时间
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t_newset;
                end
                
                %记录插入工序
                insert_index=Index1;
                                
                %更新当前机器的加工信息
                ON=ON+1;
                for j=1:ON
                    Index1=Machine(mm(Index1)).Op(j);
                    if j==1
                        Machine(mm(Index1)).GapT(j)=0;
                    else%不是第一个加工工序时
                        %对于在插入工序后的第一个工序：开始时间不变，准备时间改变-->结束时间改变，平移finish的差值
                        %对于在插入工序后的其他工序,更新工序的开始时间、完成时间：平移相同的差值
                        %若不考虑准备时间，插入到空闲位置，后面的工序开始时间不变
                        LastOp=Machine(mm(Index1)).Op(j-1);%上一个工序
                        if LastOp==insert_index%对于在插入工序后的第一个工序
                            new_finish=start{s1(Index1),s2(Index1)}+time{f_index,s1(Index1),s2(Index1),mm(Index1)};
                            shift_value=new_finish-finish{s1(Index1),s2(Index1)};
                            finish{s1(Index1),s2(Index1)}=new_finish;
                        end
                        if LastOp>insert_index%对于在插入工序后的其他工序：开始时间、完成时间：平移相同的差值
                            start{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+shift_value;
                            finish{s1(Index1),s2(Index1)}=finish{s1(Index1),s2(Index1)}+shift_value;
                        end
                        %更新空闲时间
                        Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)}-finish{s1(LastOp),s2(LastOp)};
                        if Machine(mm(Index1)).GapT(j)<0
                            Machine(mm(Index1)).GapT(j)
                        end
                    end
                    %更新每个位置的加工完成时间
                    Machine(mm(Index1)).MFT(j)=finish{s1(Index1),s2(Index1)};
                end
                
                %更新机器完工时间
                mt{mm(Index1)}=Machine(mm(Index1)).MFT(ON);
            else %√194：如果index1==0说明没有空位就要老老实实的去在后面
                %上一个加工工件：s1(Machine(mm(i)).Op(end))
                start{s1(i),s2(i)}=Machine(mm(i)).MFT(ON);%也是开始准备的时间
                %机器完工时间+准备时间
                mt{mm(i)}=start{s1(i),s2(i)}+time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(Machine(mm(i)).Op(end)),s1(i));
                finish{s1(i),s2(i)}=mt{mm(i)};
                Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
                Machine(mm(i)).GapT=[Machine(mm(i)).GapT,0];%无空闲
                Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
            end
            
        else%√203：某工件的第一个工序，且当前机器第一次加工，更新工件的开始、结束时间，机器加工的工序、完成时间
            mt{mm(i)}=time{f_index,s1(i),s2(i),mm(i)};
            start{s1(i),s2(i)}=0;
            finish{s1(i),s2(i)}=mt{mm(i)};
            Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
            Machine(mm(i)).GapT=[Machine(mm(i)).GapT,0];%空闲时间
            Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];%完成时间
        end
        
    else
        %213：开始在相同的机器上寻找合适的空位插入，如果自身的加工时间小于空闲时间则可以插入
        ON=length(Machine(mm(i)).Op);%该机器上目前的工序数
        if ON>0
            %需要修改t=time{f_index,s1(i),s2(i),mm(i)};：考虑加工时间+准备时间与空闲比较
            %t:1*ON
            t=zeros(1,ON);
            t(1)=time{f_index,s1(i),s2(i),mm(i)};%第一个空闲，作为第一个加工工序，无准备时间
            if ON>1
                for ON_index=2:ON
                    %当前加工的工序Machine(mm(i)).Op(ON_index-1)
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
                        Index1=j; %在Index工序前插入新来的工序
                        break;
                    end
                end
            end
            if Index1~=0
                %232：不是某工件的第一个工序，可插入
                Index1=Machine(mm(i)).Op(Index1);% 在s1向量中把s1(i)插入到s1(Index1)前面,Index1<i
                
                %更新s1
                tmp=s1(i);
                for j=i:-1:Index1+1
                    s1(j)=s1(j-1);
                end
                s1(Index1)=tmp;
                
                %更新s2
                tmp=s2(i);
                for j=i:-1:Index1+1
                    s2(j)=s2(j-1);
                end
                s2(Index1)=tmp;
                
                %更新加工机器mm（同上）
                tmp=mm(i);
                for j=i:-1:Index1+1
                    mm(j)=mm(j-1);
                end
                mm(Index1)=tmp;
                
                %更新当前机器的加工索引
                for j=1:ON
                    if Machine(mm(Index1)).Op(j)>=Index1
                        Machine(mm(Index1)).Op(j)=Machine(mm(Index1)).Op(j)+1;%由于前面插入工序，其他索引需要整体向后退一格
                    end
                end
                
                %更新其他所有机器的加工索引
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
                
                %更新当前机器的加工索引（待插入的工序）
                Machine(mm(Index1)).Op=[Machine(mm(Index1)).Op,Index1];
                tmp2=Machine(mm(Index1)).Op;
                Machine(mm(Index1)).Op=sort(tmp2,'ascend');
                
                %更新待插入工序的开始加工时间
                IIndex=find(Machine(mm(Index1)).Op==Index1);
                if IIndex==1%若为当前机器的第一个加工，开始时间为0
                    start{s1(Index1),s2(Index1)}=max(0,finish{s1(Index1),s2(Index1)-1});
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t(1);
                else%更新开始时间（开始准备）：max{该工件前一个工序的结束时间，当前机器前一个工序的结束时间}
                    LastOp=Machine(mm(Index1)).Op(IIndex-1);%插入的位置前一个加工工序
                    start{s1(Index1),s2(Index1)}=max(finish{s1(Index1),s2(Index1)-1},finish{s1(LastOp),s2(LastOp)});
                    %更新待插入工序的结束加工时间；考虑准备时间
                    %上一个工件：s1(LastOp)
                    t_newset=time{f_index,s1(Index1),s2(Index1),mm(Index1)}+s_time(s1(LastOp),s1(Index1));%准备时间+加工时间
                    finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t_newset;
                end
                
                
                %记录插入工序
                insert_index=Index1;
                
                %更新当前机器的加工信息
                ON=ON+1;
                for j=1:ON
                    shift_value=0;
                    Index1=Machine(mm(Index1)).Op(j);
                    if j==1
                        Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)};
                    else%不是第一个加工工序时
                        %对于在插入工序后的第一个工序：开始时间不变，准备时间改变-->结束时间改变，平移finish的差值
                        %对于在插入工序后的其他工序,更新工序的开始时间、完成时间：平移相同的差值
                        %若不考虑准备时间，插入到空闲位置，后面的工序开始时间不变
                        LastOp=Machine(mm(Index1)).Op(j-1);%上一个工序
                        if LastOp==insert_index%对于在插入工序后的第一个工序
                            new_finish=start{s1(Index1),s2(Index1)}+time{f_index,s1(Index1),s2(Index1),mm(Index1)};
                            shift_value=new_finish-finish{s1(Index1),s2(Index1)};
                            finish{s1(Index1),s2(Index1)}=new_finish;
                        end
                        if LastOp>insert_index%对于在插入工序后的其他工序：开始时间、完成时间：平移相同的差值
                            start{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+shift_value;
                            finish{s1(Index1),s2(Index1)}=finish{s1(Index1),s2(Index1)}+shift_value;
                        end
                        %更新空闲时间
                        Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)}-finish{s1(LastOp),s2(LastOp)};
                    end
                    
                    %更新每个位置的加工完成时间
                    Machine(mm(Index1)).MFT(j)=finish{s1(Index1),s2(Index1)};
                    if Machine(mm(Index1)).GapT(j)<0
                        Machine(mm(Index1)).GapT(j)
                    end
                end
                %更新机器完工时间
                mt{mm(Index1)}=Machine(mm(Index1)).MFT(ON);
            else%√290：如果index1==0说明没有空位就要老老实实的去在后面
                %不是机器上的第一个加工工序，加准备时间，上一个工件：s1(Machine(mm(i)).Op(end))
                start{s1(i),s2(i)}=max(Machine(mm(i)).MFT(ON),finish{s1(i),s2(i)-1});
                mt{mm(i)}=start{s1(i),s2(i)}+time{f_index,s1(i),s2(i),mm(i)}+s_time(s1(Machine(mm(i)).Op(end)),s1(i));
                finish{s1(i),s2(i)}=mt{mm(i)};
                Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
                gap=start{s1(i),s2(i)}-Machine(mm(i)).MFT(ON);
                Machine(mm(i)).GapT=[Machine(mm(i)).GapT,gap];
                Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
            end
        else%√300：同一个工件/当前机器第一个加工工序，不考虑准备时间
            mt{mm(i)}=finish{s1(i),s2(i)-1}+time{f_index,s1(i),s2(i),mm(i)};
            start{s1(i),s2(i)}=finish{s1(i),s2(i)-1};
            finish{s1(i),s2(i)}=mt{mm(i)};
            Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
            Machine(mm(i)).GapT=[Machine(mm(i)).GapT,start{s1(i),s2(i)}];%第一个加工工序前的空闲
            Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
        end
        
    end
end
newp=s1;
end