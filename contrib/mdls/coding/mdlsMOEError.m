function mdlsMOEError(X,D,A,lambda,b,k,mode)
    F   = X-D*A;
    F    = sum(F(:).^2);
    A    = abs(A);
    Rerr = log(A+b) - A./(A+b) + (b/(k-1))./(A+b) - (log(b)+ 1/k);
    Rerr = sum(Rerr);
    R    = sum(log(A) + b);
    switch (mode)
      case 0 % reg-constrained: nothing to do
      case 1 % error-constrained: cost function is reg alone
        T = R;
        relErr = Rerr/R;
      case 2 % lagrangian       : cost function is fit+reg
        T = F + lambda * R;
        relErr = lambda*Rerr/T;
    end
    fprintf('Absolute error=%f\trelative error=%f\n',Rerr,relErr);
end