clc
clear 
close all
DataPath='DATASET\';
job_type=[10,20,30,40,50,100,150,200];
f_type={[2],[2,3],[2,3],[2,3,4],[3,4,5],[4,5,6,7],[5,6,7],[6,7]};
nameJ='J';nameF='F'; txt='.txt';
k=1;
for i=1:size(job_type,2)
    for j=f_type{i}
        RealPath{k}=[DataPath,num2str(job_type(i)),nameJ,num2str(j),nameF,txt];
        ResultPath{k}=[num2str(job_type(i)),nameJ,num2str(j),nameF];
        k=k+1;
    end
end
global  N H SH NM M TM time AP AM AN F ProF;
File=3;
[N,F,TM,H,SH,NM,M,time,ProF]=DataReadDHFJSP(RealPath{File});

FileName=['Archive.mat'];
a = load(FileName);
AP=a.AP;
AM=a.AM;
AN=a.AN;
my_drawDFJSP(AP(1,:),AM(1,:),AN(1,:));
my_drawDFJSP(AP(2,:),AM(2,:),AN(2,:));
    