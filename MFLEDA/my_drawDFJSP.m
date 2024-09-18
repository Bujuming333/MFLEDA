%my drawDFJSP
function my_drawDFJSP(p_chrom,m_chrom,f_chrom)
global N SH F;
s1=p_chrom;
P=cell(1,F);
for i=1:SH
    t1=s1(i);%记录到当前是那个工件
    P{f_chrom(t1)+1}=[P{f_chrom(t1)+1} p_chrom(i)];
end

FJ=cell(1,F);
for i=1:N
    FJ{f_chrom(i)+1}=[FJ{f_chrom(i)+1} i];
end

%drawFJSP_SDST
for i=1:F
    drawFJSP_SDST_M5(P{i},m_chrom,FJ{i},i);
end

end