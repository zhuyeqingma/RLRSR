% 
% Usage: W=mexFistaGraph(Y,X,W0,graph,param);
%
% Name: mexFistaGraph
%
% Description: mexFistaGraph solves sparse regularized problems.
%         X is a design matrix of size m x p
%         X=[x^1,...,x^n]', where the x_i's are the rows of X
%         Y=[y^1,...,y^n] is a matrix of size m x n
%         It implements the algorithms FISTA, ISTA and subgradient descent.
%
%         It implements the algorithms FISTA, ISTA and subgradient descent for solving
%
%           min_W  loss(W) + lambda psi(W)
%          
%         The function psi are those used by mexProximalGraph (see documentation)
%         for the loss functions, see the documentation of mexFistaFlat
%         
%         This function can also handle intercepts (last row of W is not regularized),
%         and/or non-negativity constraints on W.
%
% Inputs: Y:  double dense m x n matrix
%         X:  double dense or sparse m x p matrix   
%         W0:  double dense p x n matrix or p x Nn matrix (for multi-logistic loss)
%              initial guess
%         graph: struct (see documentation of mexProximalGraph)
%         param: struct
%            param.loss (choice of loss, see above)
%            param.regul (choice of regularization, see function mexProximalFlat)
%            param.lambda (regularization parameter)
%            param.lambda2 (optional, regularization parameter, 0 by default)
%            param.lambda3 (optional, regularization parameter, 0 by default)
%            param.verbose (optional, verbosity level, false by default)
%            param.pos (optional, adds positivity constraints on the
%                coefficients, false by default)
%            param.numThreads (optional, number of threads for exploiting
%                multi-core / multi-cpus. By default, it takes the value -1,
%                which automatically selects all the available CPUs/cores).
%            param.max_it (optional, maximum number of iterations, 100 by default)
%            param.it0 (optional, frequency for computing duality gap, every 10 iterations by default)
%            param.tol (optional, tolerance for stopping criteration, which is a relative duality gap
%                if it is available, or a relative change of parameters).
%            param.gamma (optional, multiplier for increasing the parameter L in fista, 1.5 by default)
%            param.L0 (optional, initial parameter L in fista, 0.1 by default, should be small enough)
%            param.fixed_step (deactive the line search for L in fista and use param.L0 instead)
%            param.compute_gram (optional, pre-compute X^TX, false by default).
%            param.intercept (optional, do not regularize last row of W, false by default).
%            param.ista (optional, use ista instead of fista, false by default).
%            param.subgrad (optional, if not param.ista, use subradient descent instead of fista, false by default).
%            param.a, param.b (optional, if param.subgrad, the gradient step is a/(t+b)
%            also similar options as mexProximalTree
%
%            the function also implements the ADMM algorithm via an option param.admm=true. It is not documented
%            and you need to look at the source code to use it.
%
%
% Output:  W:  double dense p x n matrix or p x Nn matrix (for multi-logistic loss)
%
% Author: Julien Mairal, 2010


