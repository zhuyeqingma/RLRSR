function Y = Pre_label(gnd)
%
nClass=length(unique(gnd));
Y=zeros(nClass,length(gnd));    % 原始的标签矩阵全部为零  
for i=1:length(gnd)
    for j=1:nClass
        if j==gnd(i)
            Y(j,i)=1;   % 为有标签的样本赋标签为1
        end  
    end
end
