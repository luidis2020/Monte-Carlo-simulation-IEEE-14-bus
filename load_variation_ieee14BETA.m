if activate_normalPF==2
    
vector_PQ=ieee14_load;

end
if activate_PEM==0
vector_PQ=ieee14_load;
   
 %std14_PQ_percentage=10;
 std14_PQ=abs((std14_PQ_percentage/100).*ieee14_load);
%  weibull_k=zeros;
%  weibull_lambda=zeros;
%      for aux_column=1:2
%          for aux_row=1:length(bus)
%              weibull_k(aux_row,aux_column)=abs((std14_PQ(aux_row,aux_column)/ieee14_load(aux_row,aux_column))^(-1.086));
%              weibull_k(1,1)=0;
%              weibull_k(1,2)=0;
%              weibull_k(7,1)=0;
%              weibull_k(7,2)=0;
%              weibull_k(8,1)=0;
%              weibull_k(8,2)=0;
%              
%              weibull_lambda(aux_row,aux_column)=ieee14_load(aux_row,aux_column)/gamma(1+1/weibull_k(aux_row,aux_column));
%              pdf_PQ(aux_row,aux_column)=wblrnd(abs(weibull_lambda(aux_row,aux_column)),weibull_k(aux_row,aux_column),1);
%        %      pdf_PQ(4,2)=-1*pdf_PQ(4,2);
%          end
%      end
pdf_PQ=normrnd(ieee14_load,std14_PQ);
%pdf_PQ=ieee14_load+std14_PQ.*randn(1,1);
         
     vector_PQ(2:6,1)=pdf_PQ(2:6,1);
     vector_PQ(9:14,1)=pdf_PQ(9:14,1);
     vector_PQ(2:6,2)=pdf_PQ(2:6,2);
     vector_PQ(9:14,2)=pdf_PQ(9:14,2);
%     vector_PQ(4,2)=-1*pdf_PQ(4,2); %Weibull Distribution
     vector_PQ(4,2)=pdf_PQ(4,2);     %Normal Distribution
%vector_PQ
end