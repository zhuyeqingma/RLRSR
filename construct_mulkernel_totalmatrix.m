function K=construct_mulkernel_totalmatrix(X,Y,bet)   
% This function constructs the kernel matrix by using selected

X       =    X./ repmat(sqrt(sum(X.*X)),[size(X,1) 1]);
Y       =    Y./ repmat(sqrt(sum(Y.*Y)),[size(Y,1) 1]);

kernelparm=bet;
[d,M]=size(X);
[d,N]=size(Y);
K=zeros(M,N);
m=size(bet,1);
for i=1:1:m
    K=K+construct_mulkernel_matrix(X,Y,i) *kernelparm(i);   
end

