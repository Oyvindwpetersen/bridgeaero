function [neg_logL,neg_logL_grad]=ad_gpr_multi_loglik_wrapper(test_matrix,y,model,theta)

%%

for idx1=1:size(model,1)
    for idx2=1:size(model,2)

        theta_hyp=theta(model{idx1,idx2}.idx.glob);
        model{idx1,idx2}=theta2model(model{idx1,idx2},theta_hyp);
        
    end
end

neg_logL_grad=zeros(length(theta),1);
neg_logL=zeros(size(model,1),size(model,2));

for idx1=1:size(model,1)
    for idx2=1:size(model,2)

        [neg_logL_loc,neg_logL_grad_loc]=ad_gpr_loglik(test_matrix{idx1,idx2},y{idx1,idx2},model{idx1,idx2});
        neg_logL_grad(model{idx1,idx2}.idx.glob,1)=neg_logL_grad(model{idx1,idx2}.idx.glob,1)+neg_logL_grad_loc;
        neg_logL(idx1,idx2)=neg_logL_loc;

    end
end

neg_logL=sum(sum(neg_logL));