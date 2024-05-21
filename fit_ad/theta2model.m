function model=theta2model(model,theta)

%%

model.hyp.d=theta(model.idx.d);
model.hyp.sigma_v=theta(model.idx.sigma_v);
model.hyp.sigma=theta(model.idx.sigma);
model.hyp.L=theta(model.idx.L);
model.hyp.alpha=theta(model.idx.alpha);
model.hyp.abar=theta(model.idx.abar);

