function [f, g, H] = testfun(x)




%%

fun_obj= @(x) (x(:,1)-5).^2+x(:,2).*log(x(:,1)+15)+x(:,2).^2;


fun_grad= @(x) [ ...
                2*x(1)-10+x(2)./(x(1)+15) ;...
                log(x(1)+15)+2*x(2)];

fun_grad2= @(x) [2-x(2)./(x(1)+15).^2  1./(x(1)+15) ;
                    1./(x(1)+15)  2 ];

% Calculate objective f
f = fun_obj(x);

if nargout > 1 % gradient required
    g = fun_grad(x);

    
    if nargout > 2 % Hessian required
        H = fun_grad2(x);
    end

end