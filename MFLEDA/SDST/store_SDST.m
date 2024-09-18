%设准备时间[1,5]
job_type=[10,20,30,40,50,100,150,200];
s_time={};
for i=1:size(job_type,2)
    s_time{i}=zeros(job_type(i));
end
for i=1:size(job_type,2)
    for j=1:job_type(i)
        for k=1:job_type(i)
            if j~=k
                s_time{i}(j,k)=randperm(5,1);
            end
        end
    end
    time=s_time{i};
    file_name=['SDST_',strcat(int2str(job_type(i))),'.txt'];
    save(file_name, 'time', '-ascii');
end

% FileName1='SDST.txt';
% a=textread(FileName1);