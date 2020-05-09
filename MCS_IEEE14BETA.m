%clc
%clear all
minTime=Inf;
tic;
% 365 days per 24h
%data_amount=365*24;
%
%data_amount=5000;
%data_amount=1000;
case_IEEE14_start;
ieee14_load=barras_ipfc(:,8:9);

num_gen=0;
for num_gen_aux=1:size(barras_ipfc,1)
    if barras_ipfc(num_gen_aux,5)~=0
    num_gen=num_gen+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DN_Pcal=zeros(data_amount,size(barras_ipfc,1));
DN_Qcal=zeros(data_amount,size(barras_ipfc,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DN_P=zeros(data_amount,size(barras_ipfc,1));
DN_Q=zeros(data_amount,size(barras_ipfc,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
voltage_pdf=zeros(data_amount,size(barras_ipfc,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PgMW_pdfB=zeros(data_amount,num_gen);
QgMVAR_pdfB=zeros(data_amount,num_gen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Voltage_mean_std=zeros(size(barras_ipfc,1),2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PgQg_mean_std=zeros(num_gen,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for MCS_iter=1:data_amount
   activate_normalPF=0;
   activate_PEM=0;
   MCS_iter
   load_flow_analysis_original;
   
   DN_P(MCS_iter,:)=vector_PQ(:,1)';
   DN_Q(MCS_iter,:)=vector_PQ(:,2)';
      
   DN_Pcal(MCS_iter,:)=Pcal';
   DN_Qcal(MCS_iter,:)=Qcal';
  
%    
   PgMW_pdfB(MCS_iter,:)=gen(:,2)';
   QgMVAR_pdfB(MCS_iter,:)=gen(:,3)';
%   
  voltage_pdf(MCS_iter,:)=V';
end

for vms_aux=1:length(V)
Voltage_mean_std(vms_aux,:)=[mean(voltage_pdf(:,vms_aux)) std(voltage_pdf(:,vms_aux))];
end
for pgqg_aux=1:length(gen)
PgQg_mean_std(pgqg_aux,:)=[mean(PgMW_pdfB(:,pgqg_aux)) std(PgMW_pdfB(:,pgqg_aux)) mean(QgMVAR_pdfB(:,pgqg_aux)) std(QgMVAR_pdfB(:,pgqg_aux))];            
end
%Voltage_mean_std
tstart=tic;
telapsed=toc(tstart);
minTime=min(telapsed,minTime);
averageTimeMCS=toc;