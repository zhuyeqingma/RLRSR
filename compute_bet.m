function bet_obtain=compute_bet(Train_Ma,bet,V,Abet,Bbet,step2,R,Y,Dis,mu)
% sprintf('compute_bet_MKL')
    M=size(bet,1);
    bet_result=zeros(1,M)';

    for i=1:6
    KXX=construct_mulkernel_totalmatrix(Train_Ma,Train_Ma,bet); 
    for m=1:M
    KXXm=construct_mulkernel_matrix(Train_Ma,Train_Ma,m);
    bet_result(m) = 0.5*trace(KXXm'*V*V'*KXX+KXX'*V*V'*KXXm-2*KXXm'*V*R'*Y)+trace(V'*KXXm*Dis*KXX'*V+V'*KXX*Dis*KXXm'*V)+mu/2*trace(V'*KXXm*(Bbet*Bbet')*KXX'*V+V'*KXX*(Bbet*Bbet')*KXXm'*V-2*V'*KXXm*Bbet*Abet');
    end
    t=sqrt(sum(bet_result.^2));
    bet_result=bet+step2*bet_result/t;  
    
    row=size(bet_result,1);
    H=2*eye(row);
    f=-1*bet_result';
    Aeq=ones(1,row);
    beq= 1;
    LB=0;
    UB=ones(row,1);
    opts = optimset('Algorithm','interior-point-convex','Display','off');  %   interior-point
    bet_obtain=quadprog(H,f,[],[],Aeq,beq,LB,UB,[],opts);
    end
