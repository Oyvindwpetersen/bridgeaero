%%

nx=10
na=200


S=zeros(nx*na);
for k1=1:(nx*na)

    for k2=1:(nx*na)

        % k1_check=mod(k2+nx-1,nx)*na+ceil([k2]/nx);

        j=k2;
        q=nx;
        % k1_check=j-1-q*floor((j-1)/q)+ceil(j/q);
        k1_check=j+q-1-q*floor((j+q-1)/q)+ceil(j/q);


        if k1==k1_check
            S(k1,k2)=1;
        end
    end
end

S2=restack_a(na,nx);


dS=sparse(S-S2)

