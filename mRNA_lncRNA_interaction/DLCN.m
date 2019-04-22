clear all;
clc;
load end_allhypertest
num_mRNA=tabulate(end_allhypertest(:,1));
num_lncRNA=tabulate(end_allhypertest(:,2));
%
HSC_lnc=importdata('HSC_LNC.xlsx');
HSC_mRNA=importdata('HSC_mRNA.xlsx');
HSC_mRNA=[HSC_mRNA.textdata num2cell(HSC_mRNA.data)];
HSC_lnc=[HSC_lnc.textdata num2cell(HSC_lnc.data)];
%
num_mRNA=tabulate(end_allhypertest(:,1));
num_lncRNA=tabulate(end_allhypertest(:,2));
%
num_mRNA(:,1)= strrep(num_mRNA(:,1),'''','');
num_lncRNA(:,1)= strrep(num_lncRNA(:,1),'''','');
end_allhypertest(:,1:2)= strrep(end_allhypertest(:,1:2),'''','');
%
hsc_mRNA=[];
for i=1:size(num_mRNA,1)
    [m1,n1]=ismember(HSC_mRNA(:,1),num_mRNA{i,1}); 
     pcc_mRNA=HSC_mRNA(find(n1==1),:);
     hsc_mRNA=[hsc_mRNA;pcc_mRNA];
end
%
hsc_lncRNA=[];
for i=1:size(num_lncRNA,1)
    [m2,n2]=ismember(HSC_lnc(:,1),num_lncRNA{i,1}); 
     pcc_lncRNA=HSC_lnc(find(n2==1),:);
     hsc_lncRNA=[hsc_lncRNA;pcc_lncRNA];
end
%% 

 pcc_normal_allhypertest=[];
 for i=1:size(end_allhypertest,1)
     [m3,n3]=ismember(hsc_mRNA(:,1),end_allhypertest{i,1});
     p_mRNA=hsc_mRNA(find(n3==1),:);
     [m4,n4]=ismember(hsc_lncRNA(:,1),end_allhypertest{i,2});
     p_lncRNA=hsc_lncRNA(find(n4==1),:);
     pcc_normal=corr(cell2mat(p_mRNA(:,2:4)'),cell2mat(p_lncRNA(:,2:4)'),'type','pearson');
     pcc_normal_allhypertest=[pcc_normal_allhypertest;pcc_normal];
 end
 %% 
  pcc_cp_allhypertest=[];
 for i=1:size(end_allhypertest,1)
     [m3,n3]=ismember(hsc_mRNA(:,1),end_allhypertest{i,1});
     p_mRNA=hsc_mRNA(find(n3==1),:);
     [m4,n4]=ismember(hsc_lncRNA(:,1),end_allhypertest{i,2});
     p_lncRNA=hsc_lncRNA(find(n4==1),:);
     pcc_cp=corr(cell2mat(p_mRNA(:,5:10)'),cell2mat(p_lncRNA(:,5:10)'),'type','pearson');
     pcc_cp_allhypertest=[pcc_cp_allhypertest;pcc_cp];
 end
 %% 
   pcc_ap_allhypertest=[];
 for i=1:size(end_allhypertest,1)
     [m3,n3]=ismember(hsc_mRNA(:,1),end_allhypertest{i,1});
     p_mRNA=hsc_mRNA(find(n3==1),:);
     [m4,n4]=ismember(hsc_lncRNA(:,1),end_allhypertest{i,2});
     p_lncRNA=hsc_lncRNA(find(n4==1),:);
     pcc_ap=corr(cell2mat(p_mRNA(:,11:14)'),cell2mat(p_lncRNA(:,11:14)'),'type','pearson');
     pcc_ap_allhypertest=[pcc_ap_allhypertest;pcc_ap];
 end
 %% 
 bc=[];
    pcc_bc_allhypertest=[];
 for i=1:size(end_allhypertest,1)
     [m3,n3]=ismember(hsc_mRNA(:,1),end_allhypertest{i,1});
     p_mRNA=hsc_mRNA(find(n3==1),:);
     [m4,n4]=ismember(hsc_lncRNA(:,1),end_allhypertest{i,2});
     p_lncRNA=hsc_lncRNA(find(n4==1),:);
     bc=[bc;p_mRNA(:,15:16)' p_lncRNA(:,15:16)'];
     pcc_bc=corr(cell2mat(p_mRNA(:,15:16)'),cell2mat(p_lncRNA(:,15:16)'),'type','pearson');
     pcc_bc_allhypertest=[pcc_bc_allhypertest;pcc_bc];
 end
 %% 
 
  xlswrite('pcc_normal.xls',pcc_normal_allhypertest);
  xlswrite('pcc_cp.xls',pcc_cp_allhypertest);
  xlswrite('pcc_ap.xls',pcc_ap_allhypertest);
  xlswrite('pcc_bc.xls',pcc_bc_allhypertest);
  %% 
  normal_row=find(pcc_normal_allhypertest>0.5);
  normal_interaction_data=pcc_normal_allhypertest(normal_row,:);
  normal_interaction_gene=end_allhypertest(normal_row,1:2);
  %
  cp_row=find(pcc_cp_allhypertest>0.5);
  cp_interaction_data=pcc_cp_allhypertest(cp_row,:);
  cp_interaction_gene=end_allhypertest(cp_row,1:2);
  %
  ap_row=find(pcc_ap_allhypertest>0.5);
  ap_interaction_data=pcc_ap_allhypertest(ap_row,:);
  ap_interaction_gene=end_allhypertest(ap_row,1:2);
  %
    bc_row=find(pcc_bc_allhypertest>0.5);
  bc_interaction_data=pcc_bc_allhypertest(bc_row,:);
  bc_interaction_gene=end_allhypertest(bc_row,1:2);
  %% 
  normal_degree_m=tabulate(normal_interaction_gene(:,1));
  normal_degree_lnc=tabulate(normal_interaction_gene(:,2));
  %
  cp_degree_m=tabulate(cp_interaction_gene(:,1));
  cp_degree_lnc=tabulate(cp_interaction_gene(:,2));
  %
    ap_degree_m=tabulate(ap_interaction_gene(:,1));
  ap_degree_lnc=tabulate(ap_interaction_gene(:,2));
  %
bc_degree_m=tabulate(bc_interaction_gene(:,1));
  bc_degree_lnc=tabulate(bc_interaction_gene(:,2));
%%½»¼¯Ç°
  %% 
  xlswrite('normal_interaction.xls',[normal_interaction_gene num2cell(normal_interaction_data)]);
  xlswrite('cp_interaction.xls',[cp_interaction_gene num2cell(cp_interaction_data)]);
  xlswrite('ap_interaction.xls',[ap_interaction_gene num2cell(ap_interaction_data)]);
  xlswrite('bc_interaction.xls',[bc_interaction_gene num2cell(bc_interaction_data)]);
