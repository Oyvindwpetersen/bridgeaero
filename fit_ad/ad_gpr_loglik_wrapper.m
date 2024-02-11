function [neg_logL,neg_logL_grad]=ad_gpr_loglik_wrapper(test_matrix,y,model,theta)

%%

model.hyp.d=theta(model.idx.d);
model.hyp.sigma_v=theta(model.idx.sigma_v);
model.hyp.sigma=theta(model.idx.sigma);
model.hyp.L=theta(model.idx.L);
model.hyp.alpha=theta(model.idx.alpha);
model.abar=theta(model.idx.abar);

[neg_logL,neg_logL_grad]=ad_gpr_loglik(test_matrix,y,model);


