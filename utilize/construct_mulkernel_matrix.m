function K=construct_mulkernel_matrix(X,Y,parkm)   
% This function constructs the kernel matrix by using selected

[d,N]=size(Y);
[d,M]=size(X);

switch parkm
          %------------ 3 linear kernel -------------------
    case  1
          K=slmetric_pw(X, Y,'dotprod');   % kernel    x'x  
    case  2
          K=slmetric_pw(X, Y,'dotprod');   % kernel    x'x  
          K=0.5*ones(size(K,1),size(K,2))+K;
          % K=kernel_option.par*ones(size(K,1),size(K,2))+K;        
    case  3
          K=slmetric_pw(X, Y,'dotprod');   % kernel    x'x  
          K=1*ones(size(K,1),size(K,2))+K;
          
          %-------------- 3 polynomial kernel ------------
    case  4
          K=X'*Y;
          K=0.3*ones(size(K,1),size(K,2))+K;       %  ywz   (a+x'x)^b
          K=K.^2;
    case  5
          K=X'*Y;
          K=0.2*ones(size(K,1),size(K,2))+K;       %  ywz   (a+x'x)^b
          K=K.^3;
    case  6
          K=X'*Y;
          K=0.1*ones(size(K,1),size(K,2))+K;       %  ywz   (a+x'x)^b
          K=K.^4;
          %-------------- 6 gaussian kernel  -------------------------------
    case  7
          K=sqrDist(X,Y)';
          K=exp(-K/(2*(3^2)));
    case  8
          K=sqrDist(X,Y)';
          K=exp(-K/(2*(3.5^2)));   
    case  9
          K=sqrDist(X,Y)';
          K=exp(-K/(2*(4^2)));
          
   case  10
          K=X'*Y;
          K=0.3*ones(size(K,1),size(K,2))+K;       %  ywz   (a+x'x)^b
          K=K.^3;
    case  11
          K=X'*Y;
          K=0.4*ones(size(K,1),size(K,2))+K;       %  ywz   (a+x'x)^b
          K=K.^3;
    case  12
          K=slmetric_pw(X, Y,'dotprod');   % kernel    x'x  
          K=0.3*ones(size(K,1),size(K,2))+K;
          
          
          
          
%     case  10
%           K=sqrDist(X,Y)';
%           K=exp(-K/(2*(4.5^2)));
                   
end


