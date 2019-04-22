%Hypergeometric test screening of p>0.05 gene pairs without gene expression data
clear all;
clc;
SDE_mRNA=importdata('SDE_mRNA.xlsx');
SDE_lncRNA=importdata('SDE_lnc.xlsx');
data1=importdata('mRNA-miRNA.xlsx');  
data2=importdata('mRNA-miRNA2.xlsx');  
datalnc=importdata('lncRNA-miRNA.xlsx');  

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
%%
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
%%
interaction_mRNA=[interaction_mRNA1;interaction_mRNA3];
num_miRNAm=tabulate(interaction_mRNA(:,1));
num_mRNA=tabulate(interaction_mRNA(:,2));
num_miRNAlnc=tabulate(interaction_lncRNA(:,1));
num_lncRNA=tabulate(interaction_lncRNA(:,2));

% intersect_NUMmiRNA=intersect(num_miRNAm(:,1), num_miRNAlnc(:,1));
%num_mRNA1= strrep(num_mRNA(:,1),'''','');
union_NUMmiRNA=union(num_miRNAm(:,1), num_miRNAlnc(:,1));
Total=size(union_NUMmiRNA,1);
allhypertest=[];
for i=1:size(num_mRNA)
    for j=1:size(num_lncRNA)
        miRNA_mRNA=[];
        miRNA_lncRNA=[];
        NmRNA=num_mRNA(i,2);
        Nlnc=num_lncRNA(j,2);
        p=min(NmRNA{1,1},Nlnc{1,1});
        %Calculate the number of mRNAs and lncRNAs competing for the same miRNA
        [m1,n1]=ismember(interaction_mRNA(:,2),num_mRNA(i,1)); 
        miRNA_mRNA=interaction_mRNA(find(n1==1),:);
        %
        [m2,n2]=ismember(interaction_lncRNA(:,2),num_lncRNA(j,1)); 
        miRNA_lncRNA=interaction_lncRNA(find(n2==1),:);
        %Calculate the number of competing for the same miRNA.
        q=0;
        for ii=1:size(miRNA_mRNA,1)
            for jj=1:size(miRNA_lncRNA,1)
                a=isequal(miRNA_mRNA(ii,1),miRNA_lncRNA(jj,1));
                if a==1
                    q=q+1;
                end
            end
        end
        if q>0
            htest=[];
            for iii=q:p
                f1=combntns(NmRNA{1,1},q);
                f2=combntns(Total-NmRNA{1,1},Nlnc{1,1}-q);
                f3=combntns(Total,Nlnc{1,1});
                f=f1*f2/f3;
                htest=[htest;f];
            end
            pvalue=sum(htest);
            hypertest=[num_mRNA(i,1) num_lncRNA(j,1) pvalue];
            allhypertest=[allhypertest;hypertest];
        end
    end
end
[x,y]=find(cell2mat(allhypertest(:,3))<0.05);
end_allhypertest=(allhypertest(x,:));

%%
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
%去除元胞数组引号
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

 pcc_normal_allhypertest=[];
 for i=1:size(end_allhypertest,1)
     [m3,n3]=ismember(hsc_mRNA(:,1),end_allhypertest{i,1});
     p_mRNA=hsc_mRNA(find(n3==1),:);
     [m4,n4]=ismember(hsc_lncRNA(:,1),end_allhypertest{i,2});
     p_lncRNA=hsc_lncRNA(find(n4==1),:);
     pcc_normal=corr(cell2mat(p_mRNA(:,2:4)'),cell2mat(p_lncRNA(:,2:4)'),'type','pearson');
     pcc_normal_allhypertest=[pcc_normal_allhypertest;pcc_normal];
 end
 
  pcc_cp_allhypertest=[];
 for i=1:size(end_allhypertest,1)
     [m3,n3]=ismember(hsc_mRNA(:,1),end_allhypertest{i,1});
     p_mRNA=hsc_mRNA(find(n3==1),:);
     [m4,n4]=ismember(hsc_lncRNA(:,1),end_allhypertest{i,2});
     p_lncRNA=hsc_lncRNA(find(n4==1),:);
     pcc_cp=corr(cell2mat(p_mRNA(:,5:10)'),cell2mat(p_lncRNA(:,5:10)'),'type','pearson');
     pcc_cp_allhypertest=[pcc_cp_allhypertest;pcc_cp];
 end
 
   pcc_ap_allhypertest=[];
 for i=1:size(end_allhypertest,1)
     [m3,n3]=ismember(hsc_mRNA(:,1),end_allhypertest{i,1});
     p_mRNA=hsc_mRNA(find(n3==1),:);
     [m4,n4]=ismember(hsc_lncRNA(:,1),end_allhypertest{i,2});
     p_lncRNA=hsc_lncRNA(find(n4==1),:);
     pcc_ap=corr(cell2mat(p_mRNA(:,11:14)'),cell2mat(p_lncRNA(:,11:14)'),'type','pearson');
     pcc_ap_allhypertest=[pcc_ap_allhypertest;pcc_ap];
 end
  
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
  %
   HSC_normal_cluster=[];
  for i=1:size(normal_degree_m,1)
      [p1,q1]=ismember(hsc_mRNA(:,1),normal_degree_m{i,1});
      normal_cluster=hsc_mRNA(find(q1==1),:);
      HSC_normal_cluster=[HSC_normal_cluster;normal_cluster];
  end
    for i=1:size(normal_degree_lnc)
      [p2,q2]=ismember(hsc_lncRNA(:,1),normal_degree_lnc{i,1});
      normal_cluster=hsc_lncRNA(find(q2==1),:);
      HSC_normal_cluster=[HSC_normal_cluster;normal_cluster];
    end
    %
       HSC_cp_cluster=[];
  for i=1:size(cp_degree_m,1)
      [p1,q1]=ismember(hsc_mRNA(:,1),cp_degree_m{i,1});
      cp_cluster=hsc_mRNA(find(q1==1),:);
      HSC_cp_cluster=[HSC_cp_cluster;cp_cluster];
  end
    for i=1:size(cp_degree_lnc)
      [p2,q2]=ismember(hsc_lncRNA(:,1),cp_degree_lnc{i,1});
      cp_cluster=hsc_lncRNA(find(q2==1),:);
      HSC_cp_cluster=[HSC_cp_cluster;cp_cluster];
    end
    %
  HSC_ap_cluster=[];
  for i=1:size(ap_degree_m,1)
      [p1,q1]=ismember(hsc_mRNA(:,1),ap_degree_m{i,1});
      ap_cluster=hsc_mRNA(find(q1==1),:);
      HSC_ap_cluster=[HSC_ap_cluster;ap_cluster];
  end
    for i=1:size(ap_degree_lnc)
      [p2,q2]=ismember(hsc_lncRNA(:,1),ap_degree_lnc{i,1});
      ap_cluster=hsc_lncRNA(find(q2==1),:);
      HSC_ap_cluster=[HSC_ap_cluster;ap_cluster];
    end
    
      HSC_bc_cluster=[];
  for i=1:size(bc_degree_m,1)
      [p1,q1]=ismember(hsc_mRNA(:,1),bc_degree_m{i,1});
      bc_cluster=hsc_mRNA(find(q1==1),:);
      HSC_bc_cluster=[HSC_bc_cluster;bc_cluster];
  end
    for i=1:size(bc_degree_lnc)
      [p2,q2]=ismember(hsc_lncRNA(:,1),bc_degree_lnc{i,1});
      bc_cluster=hsc_lncRNA(find(q2==1),:);
      HSC_bc_cluster=[HSC_bc_cluster;bc_cluster];
    end
    %%
    xx=HSC_normal_cluster(:,2:size(HSC_normal_cluster,2));
%%
xx=HSC_cp_cluster(:,2:size(HSC_cp_cluster,2));
%%
xx=HSC_ap_cluster(:,2:size(HSC_ap_cluster,2));
%%
xx=HSC_bc_cluster(:,2:size(HSC_bc_cluster,2));
%%
xx=cell2mat(xx);
[number, row]=size(xx);   
 yy=pdist(xx,'correlation'); 
 zz=linkage(yy,'average'); 
 [h,t]=dendrogram(zz,10);
 cc=cophenet(zz,yy);
 %%
 Ncluster=input('Input number\n');    
c=cluster( zz,'maxclust', Ncluster );    
if(Ncluster>12)  
     Color = linspecer( Ncluster );  
     else  
     Color = linspecer( Ncluster, 'qualitative' );  
end  
 for i=1:Ncluster  
     for j = 1:number  
     if(c(j) == i)      
     hold on  
     plot(xx(j,1),xx(j,2),'o','MarkerFaceColor',Color(i,:),'MarkerEdgeColor',Color(i,:))  
     end  
     end  
 end  
 dendrogram(zz);

fid=fopen('F:\yan\GSE47927\HSC_DEG\HSC_cp_cluster_later.txt','w');  %open the txt
 
 for i = 1:Ncluster%Store the row coordinates of each category after clustering in .txt
     HCi = find(c == i);
     HCI = HCi';
     fprintf(fid,'%d ',HCI);
     fprintf(fid,'\r\n');
 end
fclose(fid);
all_data=textread('F:\yan\GSE47927\HSC_DEG\HSC_cp_cluster_later.txt');
%
aa=sum(all_data~=0,2);
cols=find(all(aa>0,2));
cluster_data=all_data(cols,:);
m=size(cluster_data,1);
%The SD of the calculated molecular concentration was significantly higher than the other groups.
fid_sdDNB=[];
fid_pccDNB=[];

fid_pccinout=fopen('pcc_inout.txt','w');  %open the txt

for ii=1:m
    Ncols=cluster_data(ii,:);
    Ncols(:,Ncols(1,:)<1)=[];

    DNB_cluster=xx(Ncols,:);
    DNB_sd=std(DNB_cluster(:));  
    fid_sdDNB=[fid_sdDNB;DNB_sd];
    
    pccall_DNB=corr(DNB_cluster');
    [g1,h1]=size(pccall_DNB);
    sumpccall_dnb=sum(sum(abs(pccall_DNB)));
    avepcc_DNB=sumpccall_dnb/(g1*h1);
    fid_pccDNB=[fid_pccDNB;avepcc_DNB];

    OTHER_cluster=xx;
    OTHER_cluster(Ncols,:)=[];
    k=size(DNB_cluster,1);
    kk=size(OTHER_cluster,1);
    %opcc
    for i=1:k                                           
        for j=1:kk
            pcc_inout=corr(DNB_cluster(i,:)',OTHER_cluster(j,:)');
            fprintf(fid_pccinout,'%d ',pcc_inout);
        end
        fprintf(fid_pccinout,'\r\n');
    end     
end
fclose(fid_pccinout);
fid_pccinout=textread('pcc_inout.txt');  %open the txt

fid_inout=[];
bb=sum(cluster_data~=0,2);
for i=1:size(bb,1)
    p=sum(bb(1:i));
    q=p-bb(i)+1;
    inout=fid_pccinout(q:p,:);
    pccinout=inout(:,any(inout));
    [k1,k2]=size(pccinout);
    sumpcc_inout=sum(sum(abs(pccinout)));
    avepcc_inout=sumpcc_inout/(k1*k2);
    fid_inout=[fid_inout;avepcc_inout];
end

SDd=fid_sdDNB;
PCCd=fid_pccDNB;
PCCo=fid_inout;
I=[];
for i=1:size(SDd,1)
    II=SDd(i)*PCCd(i)/PCCo(i);
    I=[I;II];
end
a_end=[I SDd PCCd PCCo];%small small small big
%%
%HSC_normal_DNB
DNB_matrix=[19,23,58,107,133,149,176,200,227,235,282,319,322,335,372,449,473,488,548,591,676,726,772,774,793,810,822,840,844,865,891,976,1003,1067,1081,1147,1188,1265,1270,1272,1291,1338,1397,1419,1537,1561,1594,1596,1610,1662,1676,1681,1692,1701,1716,1725,1755,1813,1822,1825,1830,1937,2005];
name=HSC_normal_cluster(:,1);
normal_normal_dnb_name=name(DNB_matrix,:);
hsc_normal_dnb_interaction=[];
normal_interaction=[normal_interaction_gene num2cell(normal_interaction_data)];
for i=size(normal_normal_dnb_name,1):size(normal_normal_dnb_name,1)
    [m1,n1]=ismember(normal_interaction(:,2),normal_normal_dnb_name{i,1}); 
     dnb_interaction=normal_interaction(find(n1==1),:);
     hsc_normal_dnb_interaction=[hsc_normal_dnb_interaction;dnb_interaction];
end
%
hsc_normal_dnb_interaction_end=[];
for i=1:size(normal_normal_dnb_name,1)-1
    [m1,n1]=ismember(hsc_normal_dnb_interaction(:,1),normal_normal_dnb_name{i,1}); 
     dnb_interaction=hsc_normal_dnb_interaction(find(n1==1),:);
     hsc_normal_dnb_interaction_end=[hsc_normal_dnb_interaction_end;dnb_interaction];
end
%
hsc_normal_dnb_degree=[];
normal_degree=[normal_degree_m;normal_degree_lnc];
for i=1:size(normal_normal_dnb_name,1)
    [m1,n1]=ismember(normal_degree(:,1),normal_normal_dnb_name{i,1}); 
     dnb_degree=normal_degree(find(n1==1),:);
     hsc_normal_dnb_degree=[hsc_normal_dnb_degree;dnb_degree];
end
xlswrite('F:\yan\GSE47927\DNB\HSC\hsc_normal_dnb_interaction_end.xlsx',hsc_normal_dnb_interaction_end);
xlswrite('F:\yan\GSE47927\DNB\HSC\hsc_normal_dnb_degree.xlsx',hsc_normal_dnb_degree);
%%
%HSC_cp_DNB
DNB_matrix=[19,43,107,120,185,201,256,380,477,506,537,605,632,661,665,682,984,1094,1115,1146,1179,1223,1249,1295,1304,1318,1390,1515];
name=HSC_cp_cluster(:,1);
hsc_cp_dnb_name=name(DNB_matrix,:);
hsc_cp_dnb_interaction=[];
cp_interaction=[cp_interaction_gene num2cell(cp_interaction_data)];
for i=size(hsc_cp_dnb_name,1)-1:size(hsc_cp_dnb_name,1)
    [m1,n1]=ismember(cp_interaction(:,2),hsc_cp_dnb_name{i,1}); 
     dnb_interaction=cp_interaction(find(n1==1),:);
     hsc_cp_dnb_interaction=[hsc_cp_dnb_interaction;dnb_interaction];
end
%
hsc_cp_dnb_interaction_end=[];
for i=1:size(hsc_cp_dnb_name,1)-2
    [m1,n1]=ismember(hsc_cp_dnb_interaction(:,1),hsc_cp_dnb_name{i,1}); 
     dnb_interaction=hsc_cp_dnb_interaction(find(n1==1),:);
     hsc_cp_dnb_interaction_end=[hsc_cp_dnb_interaction_end;dnb_interaction];
end
%
hsc_cp_dnb_degree=[];
cp_degree=[cp_degree_m;cp_degree_lnc];
for i=1:size(hsc_cp_dnb_name,1)
    [m1,n1]=ismember(cp_degree(:,1),hsc_cp_dnb_name{i,1}); 
     dnb_degree=cp_degree(find(n1==1),:);
     hsc_cp_dnb_degree=[hsc_cp_dnb_degree;dnb_degree];
end
xlswrite('F:\yan\GSE47927\DNB\HSC\hsc_cp_dnb_interaction_end.xlsx',hsc_cp_dnb_interaction_end);
xlswrite('F:\yan\GSE47927\DNB\HSC\hsc_cp_dnb_degree.xlsx',hsc_cp_dnb_degree);
%%
%HSC_ap_DNB
DNB_matrix=[13,39,40,87,116,119,142,151,167,186,190,191,215,273,326,383,419,468,472,562,651,702,733,760,782,783,800,848,852,950,1027,1101,1141,1163,1170,1193,1274,1281,1365,1379,1417,1431,1454,1497,1523,1546,1556,1600,1606,1618,1682,1697,1731,1735,1760,1774,1786,1808,1846,1859,1867,1896];
name=HSC_ap_cluster(:,1);
hsc_ap_dnb_name=name(DNB_matrix,:);
hsc_ap_dnb_interaction=[];
ap_interaction=[ap_interaction_gene num2cell(ap_interaction_data)];
for i=size(hsc_ap_dnb_name,1):size(hsc_ap_dnb_name,1)
    [m1,n1]=ismember(ap_interaction(:,2),hsc_ap_dnb_name{i,1}); 
     dnb_interaction=ap_interaction(find(n1==1),:);
     hsc_ap_dnb_interaction=[hsc_ap_dnb_interaction;dnb_interaction];
end
%
hsc_ap_dnb_interaction_end=[];
for i=1:size(hsc_ap_dnb_name,1)-1
    [m1,n1]=ismember(hsc_ap_dnb_interaction(:,1),hsc_ap_dnb_name{i,1}); 
     dnb_interaction=hsc_ap_dnb_interaction(find(n1==1),:);
     hsc_ap_dnb_interaction_end=[hsc_ap_dnb_interaction_end;dnb_interaction];
end
%
hsc_ap_dnb_degree=[];
ap_degree=[ap_degree_m;ap_degree_lnc];
for i=1:size(hsc_ap_dnb_name,1)
    [m1,n1]=ismember(ap_degree(:,1),hsc_ap_dnb_name{i,1}); 
     dnb_degree=ap_degree(find(n1==1),:);
     hsc_ap_dnb_degree=[hsc_ap_dnb_degree;dnb_degree];
end
xlswrite('F:\yan\GSE47927\DNB\HSC\hsc_ap_dnb_interaction_end.xlsx',hsc_ap_dnb_interaction_end);
xlswrite('F:\yan\GSE47927\DNB\HSC\hsc_ap_dnb_degree.xlsx',hsc_ap_dnb_degree);
%%
%HSC_BC_DNB
DNB_matrix=[25,48,124,131,132,189,257,262,295,316,357,358,377,378,537,677,723,783,787,828,867,894,956,958,1000,1021,1022,1023,1028,1046,1055,1079,1096,1182,1199,1218,1269,1331,1395,1415,1469,1476,1478,1520,1585,1603,1654,1656,1680,1740,1747,1810,1830,1836,1843,1882,1919,1927,2066,2105,2110,2117,2199,2282,2378,2380,2536,2567,2587,2601,2603];
name=HSC_bc_cluster(:,1);
normal_bc_dnb_name=name(DNB_matrix,:);
hsc_bc_dnb_interaction=[];
bc_interaction=[bc_interaction_gene num2cell(bc_interaction_data)];
for i=size(normal_bc_dnb_name,1):size(normal_bc_dnb_name,1)
    [m1,n1]=ismember(bc_interaction(:,2),normal_bc_dnb_name{i,1}); 
     dnb_interaction=bc_interaction(find(n1==1),:);
     hsc_bc_dnb_interaction=[hsc_bc_dnb_interaction;dnb_interaction];
end
%
hsc_bc_dnb_interaction_end=[];
for i=1:size(normal_bc_dnb_name,1)-1
    [m1,n1]=ismember(hsc_bc_dnb_interaction(:,1),normal_bc_dnb_name{i,1}); 
     dnb_interaction=hsc_bc_dnb_interaction(find(n1==1),:);
     hsc_bc_dnb_interaction_end=[hsc_bc_dnb_interaction_end;dnb_interaction];
end
%
hsc_bc_dnb_degree=[];
bc_degree=[bc_degree_m;bc_degree_lnc];
for i=1:size(normal_bc_dnb_name,1)
    [m1,n1]=ismember(bc_degree(:,1),normal_bc_dnb_name{i,1}); 
     dnb_degree=bc_degree(find(n1==1),:);
     hsc_bc_dnb_degree=[hsc_bc_dnb_degree;dnb_degree];
end
xlswrite('F:\yan\GSE47927\DNB\HSC\hsc_bc_dnb_interaction_end.xlsx',hsc_bc_dnb_interaction_end);
xlswrite('F:\yan\GSE47927\DNB\HSC\hsc_bc_dnb_degree.xlsx',hsc_bc_dnb_degree);