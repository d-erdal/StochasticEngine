

% For some reason, when starting matlab new, STK needs to be run with the
% traning we plan for it to evluate the model correctly. Hence, this is
% just a dummy code

R=rand(100,32);
obs=randn(100,1);

% dimension
dim = size(R,2);

% Parameters are preset in the proxySetup-stucture
SIGMA2 = 1.0;  % variance parameter
NU     = 2.0;  % regularity parameter
RHO1   = 0.4*ones(dim,1);  % scale (range) parameter
param0 = log ([SIGMA2; 1./RHO1]);

% initialte the GPE-model
GPE_XX.model = stk_model ('stk_materncov32_aniso', dim);

% model.lm = stk_lm_constant;
GPE_XX.model.lm = stk_lm_affine;

GPE_XX.model.param = stk_param_estim (GPE_XX.model, R, obs, param0);

% Setup the model
GPE_XX.M_post = stk_model_gpposterior(GPE_XX.model, R, obs);

stk_predict (GPE_XX.M_post ,rand(1,dim));

clear GPE_XX R obs param0 RHO1 NU SIGMA2 dim










