%miRNA_mRNA_union.m is the union of miRNA_mRNA1 and miRNA_mRNA2
clear all;
clc;
interaction_mRNA1=importdata('interaction_mRNA1.xlsx');
interaction_mRNA2=importdata('interaction_mRNA2.xlsx');
m1=size(interaction_mRNA1,1);
interaction_mRNA3=interaction_mRNA2;
for i=1:m1
    for j=1:size(interaction_mRNA3,1)
        a=isequal(interaction_mRNA1(i,:),interaction_mRNA3(j,:));
        if a==1
            interaction_mRNA3(j,:)=[];
            break;
        end
    end
end
miRNA_mRNA_union=[interaction_mRNA1;interaction_mRNA3];