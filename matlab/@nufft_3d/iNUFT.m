%% inverse non-uniform FT (cartesian image <- irregular kspace)
function im = iNUFT(obj,raw,maxit,damp,W,constraint,lambda)
%im = iNUFT(obj,raw,maxit,damp,W,constraint,lambda)
%
% -raw: complex raw kspace data [nr nc] or [nr ny nc]
% -maxit [scalar]: no. iterations (use 0 or 1 for regridding)
% -damp [scalar]: Tikhonov regularization (only when maxit>1)
% -W [scalar, ny or nr*ny]: weighting (only when maxit>1)
%
% Experimental options
% -constraint: 'phase-constraint'
%              'compressed-sensing'
%              'parallel-imaging-sake'
%              'parallel-imaging-pruno'
% -lambda [scalar]: regularization parameter (e.g. 0.5, 1e-3, 0, 0.5)
%
% The constraints can often produce slightly better images but require
% tweaking of lambda. They are implemented as mutually exclusive options
% for evaluation but actually they could be used simultaneously at great
% computational cost and probably only marginal benefit.
%
%% argument checks

% expected no. data points
nrow = size(obj.H,1);

if size(raw,1)==nrow
    nc = size(raw,2);
    fprintf('  %s received raw data: nr=%i nc=%i\n',mfilename,nrow,nc);
else
    nr = size(raw,1); % assume readout points
    ny = size(raw,2); % assume radial spokes
    if nr*ny ~= nrow
        error('raw data leading dimension(s) must be length %i (not %ix%i)',nrow,nr,ny)
    end
    nc = size(raw,3);
    fprintf('  %s received raw data: nr=%i ny=%i nc=%i\n',mfilename,nr,ny,nc);
end
raw = reshape(raw,nrow,nc);

% optional argument checks
if ~exist('maxit','var') || isempty(maxit)
    maxit = 10; % default to iterative
else
    validateattributes(maxit,{'numeric'},{'scalar','finite','integer','nonnegative'},'','maxit');
end
if ~exist('damp','var') || isempty(damp)
    damp = 0;
else
    validateattributes(damp,{'numeric'},{'scalar','finite','nonnegative'},'','damp');
end
if ~exist('W','var') || isscalar(W) || isempty(W)
    W = 1;
else
    if numel(W)~=nrow
        if ~exist('ny','var')
            % guess - expect W to be vector of length ny
            ny = numel(W);
            nr = nrow/ny;
        end
        % this should catch size mismatches
        if mod(nr,1) || nr*ny~=nrow
            error('W must be a vector of length ny or ny*nr');
        end
        W = repmat(reshape(W,1,ny),nr,1);
    end
    W = reshape(W,nrow,1);
    if numel(unique(W))==1; W = W(1); end
    if ~any(W); error('W cannot be all zero'); end
    validateattributes(W,{'numeric','gpuArray'},{'finite','nonnegative'},'','W');
end
if ~exist('constraint','var') || isempty(constraint)
    constraint = '';
else
    switch constraint
        case 'phase-constraint';
        case 'compressed-sensing';
        case 'parallel-imaging-sake'; if nc==1; error('sake-low-rank requires multiple coils'); end    
        case 'parallel-imaging-pruno'; if nc==1; error('parallel-imaging requires multiple coils'); end            
        otherwise; error('unknown constraint');
    end
    if ~exist('lambda','var') || isempty(lambda)
        error('lambda must be supplied with %s',constraint);
    end
    validateattributes(lambda,{'numeric'},{'scalar','finite','nonnegative'},'','lambda');
end

% damp, weighting and constraints require iterative recon
if damp~=0 && maxit<=1
    warning('damp is only effective when maxit>1 - try 10');
end
if ~isscalar(W) && maxit<=1
    warning('weighting is only effective when maxit>1 - try 10');
end
if ~isempty(constraint) && lambda && maxit<=1
    warning('phase constraint is only effective when maxit>1 - try 10');
end

%% finalize setup
fprintf('  maxit=%i damp=%.3f weighted=%i',maxit,damp,~isscalar(W));
if isempty(constraint)
    fprintf('\n')
else
    fprintf(' (%s lambda=%.3f)\n',constraint,lambda);
end

% send to gpu if needed
if obj.gpu
    W = gpuArray(W);
    raw = gpuArray(single(raw));
else
    raw = double(gather(raw));
end

%% iNUFT reconstruction: solve Ax=b
tic;

if maxit==0
    
    % regridding x = A'Db
    x = obj.aNUFT(obj.d.*raw);
    
else
    
    % linear operator (A'WDA)
    A = @(x)obj.iprojection(x,damp,W);
    
    % rhs vector b = (A'WDb)
    b = obj.aNUFT((W.*obj.d).*raw);

    % correct shape for solver
    b = reshape(b,prod(obj.N),nc);
    
    % solve (A'WDA)(x) = (A'WDb) + penalty on ||x||
    [x,~,relres,iter] = pcgpc(A,b,[],maxit);
    
end

%% experimental options

if ~isempty(constraint)
    
    % phase constrained least squares
    if isequal(constraint,'phase-constraint')
        
        % smoothing kernel (in image space)
        h = exp(-(-obj.low:obj.low).^2/obj.low);
        
        % use regridding solution for phase
        P = reshape(x,obj.N(1),obj.N(2),obj.N(3),nc);
        P = circshift(P,fix(obj.N/2)); % mitigate edge effects (or use cconvn)
        P = convn(P,reshape(h,numel(h),1,1),'same');
        P = convn(P,reshape(h,1,numel(h),1),'same');
        P = convn(P,reshape(h,1,1,numel(h)),'same');
        P = circshift(P,-fix(obj.N/2)); % shift back again
        P = exp(i*angle(P));
        
        % linear operator (P'A'WDAP)
        A = @(x)obj.pprojection(x,damp,lambda,W,P);
        
        % rhs vector b = (P'A'WDb)
        b = conj(P).*obj.aNUFT((W.*obj.d).*raw);
        
        % correct shape for solver
        P = reshape(P,prod(obj.N),nc);
        b = reshape(b,prod(obj.N),nc);
        
        % solve (P'A'WDAP)(P'x) = (P'A'WDb) + penalty on imag(P'x)
        [x,~,relres,iter] = pcgpc(A,b,[],maxit);
        
        % put back the low resolution phase
        x = P.*x;
        
    end
    
    % compressed sensing (wavelet)
    if isequal(constraint,'compressed-sensing')
        
        % wrapper to dwt/idwt (any orthogonal choice)
        Q = DWT(obj.N,'db2'); % Q=forward Q'=inverse
        
        % rhs vector b = (QA'WDb)
        b = Q*obj.aNUFT((W.*obj.d).*raw);
        
        % correct shape for solver
        b = reshape(b,[],1);
        
        % linear operator (QA'WDAQ')
        A = @(q)reshape(Q*obj.iprojection(Q'*q,damp,W),[],1);
        
        % solve (QA'WDAQ')(q) = (QA'WDb) + penalty on ||q||_1
        q = pcgL1(A,b,lambda);
        
        % q in wavelet domain so x = Q' * q
        x = Q' * q;
        
    end
    
    % parallel imaging (sake low rank)
    if isequal(constraint,'parallel-imaging-sake')
        
        % loraks is either on or off
        x = obj.sake(raw,'damp',damp,'W',W,'loraks',lambda>0);
        
    end
    
    % parallel imaging (pruno low rank)
    if isequal(constraint,'parallel-imaging-pruno')
        
        % use regridding solution for calibration
        x = reshape(x,obj.N(1),obj.N(2),obj.N(3),nc);
        x = fft(fft(fft(x,[],1),[],2),[],3); % kspace
        
        % make the nulling kernels
        obj.fnull = obj.pruno(x) * lambda;
        obj.anull = conj(obj.fnull);
        
        % correct RHS for solver
        b = obj.aNUFT((W.*obj.d).*raw);
        b = reshape(b,[],1);
        
        % in Cartesian the diagonal of iprojection is the mean of
        % the sample weighting squared (used as a preconditioner)
        D = obj.dwsd + damp^2 + real(sum(obj.anull.*obj.fnull,5));
        M = @(x) x./reshape(D,size(x)); % diagonal preconditioner -
        
        % check: measure diagonal of iprojection (V. SLOW)
        if 0
            N = 200; % how many to test
            d = zeros(1,N,'like',D);
            for j = 1:N
                tmp = zeros(size(D));
                tmp(j) = 1; tmp = obj.iprojection(tmp,damp,W);
                d(j) = real(tmp(j)); fprintf('%i/%i\n',j,N);
            end
            plot([d;D(1:N);d-D(1:N)]'); legend({'exact','estimate','diff'});
            keyboard
        end
        
        % least squares (A'WA)(x) = (A'Wb) + penalty on ||null*x||
        iters = 100; % need about 100
        x = pcgpc(@(x)obj.iprojection(x,damp,W),b,[],iters,M);
        
    end
    
end

%% reshape into image format
im = reshape(gather(x),obj.N(1),obj.N(2),obj.N(3),nc);

fprintf('  %s returned %ix%ix%ix%i dataset. ',mfilename,obj.N(1),obj.N(2),obj.N(3),nc); toc;

end

function [x,flag,relres,iter,resvec] = pcgpc(A,b,tol,maxit,M1,M2,x0)
% Preconditioned conjugate gradient (pcg) with modifications to
% support the use of a penalty on Im(x) aka a phase constraint.
%
% Meant to be used with anonymous functions, A = @(x)myfunc(x),
% where myfunc(x) should return:
% - A*x           : to minimize ||b-Ax||^2
% - A*x+λ*x       : to minimize ||b-Ax||^2 + λ^2||x||^2
% - A*x+λ*i*Im(x) : to minimize ||b-Ax||^2 + λ^2||Im(x)||^2
%
% References:
% - An Introduction to the CG Method Without the Agonizing Pain
%    Jonathan Richard Shewchuk (1994)
% - Partial Fourier Partially Parallel Imaging
%    Mark Bydder and Matthew D. Robson (MRM 2005;53:1393)
%
% Modifications:
% - uses the real part only of the dot products
% - allows multiple RHS vectors (b = [b1 b2 ...])
% - mostly compatible with pcg (except M2 and flag)
%
% Usage: see Matlab's pcg function (help pcg)

% check arguments
if nargin<2; error('Not enough input arguments.'); end
if ~exist('tol') || isempty(tol); tol = 1e-6; end
if ~exist('maxit') || isempty(maxit); maxit = 20; end
if ~exist('M1') || isempty(M1); M1 = @(arg) arg; end
if exist('M2') && ~isempty(M2); error('M2 argument not supported'); end
if ~ismatrix(b); error('b argument must be a column vector or 2d array'); end
validateattributes(tol,{'numeric'},{'scalar','nonnegative','finite'},'','tol');
validateattributes(maxit,{'numeric'},{'scalar','nonnegative','integer'},'','maxit');

% not intended for matrix inputs but they can be supported
if isnumeric(A); A = @(arg) A * arg; end
if isnumeric(M1); M1 = @(arg) M1 \ arg; end

% initialize
t = tic;
iter = 1;
flag = 1;
if ~exist('x0') || isempty(x0)
    r = b;
    x = zeros(size(b),'like',b);
else
    if ~isequal(size(x0),size(b))
        error('x0 must be a column vector of length %i to match the problem size.',numel(b));
    end
    x = x0;
    r = A(x);   
    if ~isequal(size(r),size(b))
        error('A(x) must return a column vector of length %i to match the problem size.',numel(b));
    end  
    r = b - r;
end
d = M1(r);
if ~isequal(size(d),size(b))
    error('M1(x) must return a column vector of length %i to match the problem size.',numel(b));
end
delta0 = vecnorm(b);
delta_new = real(dot(r,d));
resvec(iter,:) = vecnorm(r);
solvec(iter,:) = vecnorm(x);

% min norm solution
xmin = x;
imin = zeros(1,size(b,2));

% main loop
while maxit
    
    iter = iter+1;
	clear q; q = A(d);
	alpha = delta_new./real(dot(d,q));
    
    % unsuccessful termination
    if ~all(isfinite(alpha)); flag = 4; break; end
    
	x = x + alpha.*d;
    
    % recalculate residual occasionally
    if mod(iter,20)==0
        r = b - A(x);
    else
        r = r - alpha.*q;
    end

    % residual and solution vectors
    resvec(iter,:) = vecnorm(r);
    solvec(iter,:) = vecnorm(x);

    % keep best solution
    ok = resvec(iter,:) < min(resvec(1:iter-1,:));
    if any(ok)
        xmin(:,ok) = x(:,ok);
        imin(ok) = iter;
    end
    
    % successful termination
    if all(resvec(iter,:)<tol*delta0); flag = 0; break; end

    % unsuccessful termination
    if iter>maxit; flag = 1; break; end

    clear q; q = M1(r);
	delta_old = delta_new;
    delta_new = real(dot(r,q));
   
    % unsuccessful termination
    if all(delta_new<=0); flag = 4; break; end

    beta = delta_new./delta_old;
    d = q + beta.*d;

end

% min norm solution
ok = imin==iter;
if any(~ok)
    flag = 3;
    x(:,~ok) = xmin(:,~ok);
end

% remeasure final residual
if nargout>2
    resvec(end,:) = vecnorm(b-A(x));
end
relres = resvec(end,:)./delta0;

% only display if flag not supplied
if nargout<2
    for k = 1:size(b,2)
        fprintf('pcg terminated at iteration %i (flag %i): relres = %e.\n',imin(k),flag,relres(k)); 
    end
    toc(t);
end
end
