
function [Train_Ma,Train_Lab,Test_Ma,Test_Lab,Save_Train_Ma] = gen_data_random(name, sele_num)
%name = 'YaleB_32x32'
load (name);

nnClass = length(unique(gnd));
num_Class=[];
for i = 1:nnClass
    num_Class = [num_Class length(find(gnd==i))];  
end
Train_Ma  = [];
Train_Lab = [];
Test_Ma   = [];
Test_Lab  = [];
for j = 1:nnClass
    idx=find(gnd==j); 
    randIdx  = randperm(num_Class(j));   
    Train_Ma = [Train_Ma; fea(idx(randIdx(1:sele_num)),:)];  
    Train_Lab= [Train_Lab;gnd(idx(randIdx(1:sele_num)))];
    Test_Ma  = [Test_Ma;fea(idx(randIdx(sele_num+1:num_Class(j))),:)];
    Test_Lab = [Test_Lab;gnd(idx(randIdx(sele_num+1:num_Class(j))))];
end
Train_Ma = Train_Ma';
Save_Train_Ma =Train_Ma;

Train_Ma = Train_Ma./repmat(sqrt(sum(Train_Ma.^2)),[size(Train_Ma,1) 1]);
Test_Ma  = Test_Ma';
Test_Ma  = Test_Ma./repmat(sqrt(sum(Test_Ma.^2)),[size(Test_Ma,1) 1]);





