warning('off');
clc
clear all;
clear memory;
addpath('./lib');
addpath('./utilize');
addpath('distance');
addpath('dataset');
load ORL_32x32
result=[]; 
for  result_iter=1:10
result_iter

[Train_Ma,Train_Lab,Test_Ma,Test_Lab,Save_Train_Ma]=gen_data_random('ORL_32x32',5);
% [Train_Ma,Train_Lab,Test_Ma,Test_Lab,Save_Train_Ma]=gen_data_random('COIL20',10);


lambda1 = 0.1;  
lambda2 = 0.01;  
lambda3 = 0.01; 
lambda5 =0.001;   
dim = 100;      

[betloop,V,Z,E,R,obj] = MK_LR_multiple_kernel(Train_Ma,Train_Lab,Test_Ma,lambda1,lambda2,lambda3,lambda5, dim);

Train_Ma_K=construct_mulkernel_totalmatrix(Train_Ma,Train_Ma,betloop);  
Test_Ma_K=construct_mulkernel_totalmatrix(Train_Ma,Test_Ma,betloop); 


%testing
Test_Maa  = V'*Test_Ma_K;
Test_Maa  = Test_Maa./repmat(sqrt(sum(Test_Maa.^2)),[size(Test_Maa,1) 1]);
Train_Maa = V'*Train_Ma_K;
Train_Maa = Train_Maa./repmat(sqrt(sum(Train_Maa.^2)),[size(Train_Maa,1) 1]);

mdl = fitcknn(real(Train_Maa'),Train_Lab,'NumNeighbors',2);
[pred] = predict(mdl,real(Test_Maa'));
acc_test = sum(Test_Lab == pred)/length(Test_Lab)*100;
result=[result acc_test]

end

t1=mean(result)


