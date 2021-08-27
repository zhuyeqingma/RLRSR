function [bet_obtain,V,Z,E,R,obj] = MK_LR_multiple_kernel(Train_Ma,Train_Lab,Test_Ma,lambda1,lambda2,lambda3,lambda5, dim)
options = [];
options.betax = 0.01; %lamda1_x 
options.intraKx = 1; % k1 
options.interKx = 5; % k2
options.nIter = 3; % iteration num.
options.lamda2 = 0.1; % lamda2
kernelnumber=10;   %
my_mu = 0.01;  % 0.01

betloop= 1/sqrt(kernelnumber)*ones(kernelnumber,1);

     %  Traing data kernelised
     Train_Ma_K=construct_mulkernel_totalmatrix(Train_Ma,Train_Ma,betloop);  
     KDD = Train_Ma_K;
     
     Lx =Train_Lab;
     Xx = KDD;


% sample sizes and dimensions
[dx, Nx] = size(Xx);
betax = 0.2;  %0.2, 0.10
% options.intraK = options.intraKx; 
% options.interK = options.interKx; 
[LWx, LBx] = CalAffinityMatrix(Lx, options, Xx');%Ax, Xx
LWx = (betax/Nx)*LWx;
LBx = (betax/Nx)*LBx;
for i=1:size(LBx,1)
    LBx(i,i)=0;
end
for i=1:size(LWx,1)
    LWx(i,i)=1000000;
end
DCol = full(sum(LWx,2));
D_LWx = spdiags(DCol,0,speye(size(LWx,1)));
L_LWx = D_LWx - LWx;
clear DCol
DCol = full(sum(LBx,2));
D_LBx = spdiags(DCol,0,speye(size(LBx,1)));
L_LWx = D_LBx - LBx;
Dis=Xx*(L_LWx-my_mu*L_LWx)*Xx';   % 




X=Train_Ma_K;
X_label=Train_Lab;
MYDis=Dis;


Max_iter =10;
[m,n] = size(X);
c = length(unique(X_label));
Y = Pre_label(X_label);

V = zeros(m,dim);
Z = zeros(n,n);
E = zeros(dim,n);
R = zeros(c,dim);
A = zeros(n,n);

Y1 = zeros(dim,n);
Y2 = zeros(n,n);
obj=zeros(Max_iter,1);
mu = 0.1;
rho = 1.01;
max_mu = 10^6;




for iter = 1:Max_iter
    %      iter
X=Train_Ma_K;
X_label=Train_Lab;
MYDis=Dis;


        [U1,~,S1] = svd(Y*X'*V,'econ');
        R = U1*S1';
        clear M;
    
    % V
    v  = sqrt(sum(V.*V,2)+eps);
    D  = diag(1./(v));
    M  = E-Y1/mu;
    V1 = X*X'+2*lambda1*MYDis +mu*(X-X*Z)*(X-X*Z)'+lambda5*D;

        
    V2 = mu*(X-X*Z)*M'+X*Y'*R;
    V=V1\V2;
        clear M;clear D;%clear v;

    % Z
    M = V'*X-E+Y1/mu;
    N = A-Y2/mu;
    Z1 = X'*V*V'*X+ eye(n);
    Z2 = X'*V*M+N;
    Z=Z1\Z2;
    clear M; clear N;
    
    % E
    inmu = lambda3/mu;
    temp_E = V'*X-V'*X*Z+Y1/mu;
    for i = 1:dim
        w = temp_E(i,:);
        la = sqrt(w*w');
        lam = 0;
        if la > inmu
            lam = 1-inmu/la;
        elseif la <= -inmu
            lam = 1+inmu/la;
        end;
        E(i,:) = lam*w;
    end;
    clear temp_E;
    
    %  A
    eps1 = lambda2/mu;
    temp_A = Z+Y2/mu;
    [uu,ss,vv] = svd(temp_A,'econ');
    ss = diag(ss);
    SVP = length(find(ss>eps1));
    if SVP>=1
        ss = ss(1:SVP)-eps1;
    else
        SVP = 1;
        ss = 0;
    end
    H = uu(:,1:SVP)*diag(ss)*vv(:,1:SVP)';
    clear temp_A; clear eps1;
    

    % Y1;Y2;mu
    Y1 = Y1+mu*(V'*X-V'*X*Z-E);
    Y2 = Y2+mu*(Z-A);
    mu = min(rho*mu,max_mu);
    leq1 = norm(V'*X-V'*X*Z-E,Inf);
    leq2 = norm(Z-A,Inf);
    ee  = sqrt(sum(E.*E,2));
    
    
    % bet
    Abet=E-Y1/mu;
    Bbet=eye(size(Z,1))-Z;
    step2=0.01;
    bet_obtain=compute_bet(Train_Ma,betloop,V,Abet,Bbet,step2,R,Y,Dis,mu);

    obj(iter) = norm(Y-R*V'*X,'fro')^2+lambda5*sum(v)+lambda3*sum(ee)+lambda2*rank(Z)+trace(V'*MYDis*V);
    
%    bet_obtain = betloop;
   
   
    betloop= bet_obtain;
   Train_Ma_K=construct_mulkernel_totalmatrix(Train_Ma,Train_Ma,betloop);  
   KDD = Train_Ma_K;
   Xx = KDD;
   % sample sizes and dimensions
%     [dx, Nx] = size(Xx);
%     betax = 0.2;  %0.2, 0.10
%     options.intraK = options.intraKx; 
%     options.interK = options.interKx; 
    [LWx, LBx] = CalAffinityMatrix(Lx, options, Xx');%Ax, Xx
    LWx = (betax/Nx)*LWx;
    LBx = (betax/Nx)*LBx;
    for i=1:size(LBx,1)
        LBx(i,i)=0;
    end
    for i=1:size(LWx,1)
        LWx(i,i)=1000000;
    end
    DCol = full(sum(LWx,2));
    D_LWx = spdiags(DCol,0,speye(size(LWx,1)));
    L_LWx = D_LWx - LWx;
    clear DCol
    DCol = full(sum(LBx,2));
    D_LBx = spdiags(DCol,0,speye(size(LBx,1)));
    L_LWx = D_LBx - LBx;
    Dis=Xx*(L_LWx-my_mu*L_LWx)*Xx';   % 
   

    if iter > 2
        if leq1 < 10^-5 && leq2 < 10^-5 && abs(obj(iter)-obj(iter-1)) < 10^-3
            break
        end
    end
end


