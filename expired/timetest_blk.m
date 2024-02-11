

%%
clc

clear a
for k=1:20
    a{k}=randn(30);
    % a{k}=sparse(a{k});
end

for k=1:20
    a{k}=randn(30);
end

tic
for k=1:1000
A0=blkdiag(a{:});
end
toc

tic
for k=1:1000
A1=blkdiag2(a{:});
end
toc

tic
for k=1:1000
A2=blkdiag2(a{:},'s');
end
toc
