clear all;
clc;
SDE_mRNA=importdata('SDE_mRNA.xlsx');
SDE_lncRNA=importdata('SDE_lnc.xlsx');
data1=importdata('mRNA-miRNA.xlsx');  %加载数据
data2=importdata('mRNA-miRNA2.xlsx');  %加载数据
datalnc=importdata('lncRNA-miRNA.xlsx');  %加载数据

DEG_mRNA=union(SDE_mRNA.CP,SDE_mRNA.AP);
DEG_mRNA=union(DEG_mRNA,SDE_mRNA.BC);
DEG_lncRNA=union(SDE_lncRNA.cp,SDE_lncRNA.ap);
DEG_lncRNA=union(DEG_lncRNA,SDE_lncRNA.bc);

m=size(DEG_mRNA,1);
m2=size(DEG_lncRNA,1);
p1=size(data1,1);
p2=size(data2,1);
plnc=size(datalnc,1);
interaction_mRNA1=[];
for i=1:m
    for j=2:p1
        if strcmp(DEG_mRNA{i,1},data1{j,2})==1
            interaction_mRNA1=[interaction_mRNA1;data1(j,:)];
        end
    end
end
% xlsxwrite('interaction_mRNA1',interaction_mRNA);

interaction_mRNA2=[];
for i=1:m
    for j=2:p2
        if strcmp(DEG_mRNA{i,1},data2{j,2})==1
            interaction_mRNA2=[interaction_mRNA2;data2(j,:)];
        end 
    end
end

interaction_lncRNA=[];
for i=1:m2
    for j=2:plnc
        if strcmp(DEG_lncRNA{i,1},datalnc{j,2})==1
            interaction_lncRNA=[interaction_lncRNA;datalnc(j,:)];
        end 
    end
end