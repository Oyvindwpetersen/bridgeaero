%function [Ka_blk,Ka_grad_tmp]=kernel_cov(x1,x2,model)

%%
%

%%

if nargout>1
    return_grad=true;
else
    return_grad=false;
end

Ka_blk=cell(1,model.na);
Ka_grad_tmp=cell(1,model.na);

for k=1:model.na

    if strcmpi(model.kernel,'se')

        if return_grad
            [Ka_blk{k},Ka_grad_tmp{k}]=kernel_se(x1,x2,model.sigma(k),model.L(k));
        else
            [Ka_blk{k}]=kernel_se(x1,x2,model.sigma(k),model.L(k));
        end

    elseif strcmpi(model.kernel,'matern52')

        if return_grad
            [Ka_blk{k},Ka_grad_tmp{k}]=kernel_matern52(x1,x2,model.sigma(k),model.L(k));
        else
            [Ka_blk{k}]=kernel_matern52(x1,x2,model.sigma(k),model.L(k));
        end


    elseif strcmpi(model.kernel,'rq')

        if return_grad
            [Ka_blk{k},Ka_grad_tmp{k}]=kernel_rq(x1,x2,model.sigma(k),model.L(k),model.alpha(k));
        else
            [Ka_blk{k}]=kernel_rq(x1,x2,model.sigma(k),model.L(k),model.alpha(k));
        end

    else
        error(['Kernel type not recognized:' model.kernel]);
    end

end

