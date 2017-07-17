function R = mdlsModelEnergy(X,D,A,params)
  if ~exist('params','var')
      params = mdlsDefaultModelParams();
  end
  %
  % the model energy is the corresponding objective function
  % that was minimized during sparse coding.
  % 
  % for reg_mode 0, the constraint was put on the regularizer,
  % so the energy is the reconstruction error (see scLasso)
  %
  % for reg_mode 1, the constraint is on the fitting term,
  % so the energy is the regularizer evalueted at the coefficients
  %
  % for reg_mode 2, the problem is unconstrained and ent obj
  % function is the Lagrangian fit + lambda * reg
  %
  switch (params.reg_mode)
    case 1 % energy = regularizer
      switch (params.reg_type)
        case 'l1'
          R = full( sum(abs(A)) );
        case 'wl1'
          R = full( params.theta_row'*abs(A) );
        case 'wl2'
          S = 1./(mean(A.^2,2)+eps);
          R = full( S'*abs(A) );          
        case 'moe'
          R = sum( log( abs(A) + params.beta) );
        case 'joe'
          a = params.lambda_min;
          b = params.lambda_max;
          C = -log( (b-a) / log(b/a) );
          inz = find(A);
          R = C*ones(size(A));
          R(inz) = log(abs(Anz)) - log(exp(-a*Anz) - exp(-b*Anz));
          R = E + params.lambda*sum(R);
        otherwise
          error([params.reg_type ' is not a valid reg type.']);
      end
    case 0 % energy = fitting term
      E = single(X - D*A);
      R = 0.5*sum(E.^2); 
    case 2 % energy = Lagrangian
      E = single(X - double(D)*A); 
      E = 0.5*sum(E.^2);
      switch (params.reg_type)
        case 'l2'
          R = E + params.lambda*full( sum(A.^2) );
        case 'l1'
          R = E + params.lambda*full( sum(abs(A)) );
        case 'wl1'
          R = E + params.lambda/mean(params.theta_row)*full( params.theta_row'*abs(A) );
        case 'wl2'
          S = 1./(mean(A.^2,2)+eps);
          R = E + params.lambda*full(S'*abs(A) );
        case 'moe'
          R = E + params.lambda*sum( log( abs(A) + params.beta) );
        case 'joe'
          a = params.lambda_min;
          b = params.lambda_max;
          C = -log( (b-a) / log(b/a) );
          inz = find(A);
          Anz = A(inz);
          R = C*ones(size(A));
          R(inz) = log(abs(Anz)) - log(exp(-a*Anz) - exp(-b*Anz));
          R = E + params.lambda*sum(R);
        otherwise
          error([params.reg_type ' is not a valid reg type.']);
      end % switch no reg type
      otherwise
        error([params.reg_mode ' is not a valid reg mode.']);
  end % switch on reg mode
end
