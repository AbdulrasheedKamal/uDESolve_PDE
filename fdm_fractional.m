function sol = fdm_fractional(pde_info, params)
% FDM solver for time-fractional PDE: D^alpha u = c*u_xx + f, alpha in (0,1).
% Ref: Liu, Anh, Turner (2004) ANZIAM J., Ehighie (2024)

    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core'));

    grid  = grid_generator(params);
    c     = pde_info.coeff.c;
    alpha = pde_info.coeff.alpha_frac;
    dt    = params.dt;
    x     = grid.x;  Nx=grid.Nx;  dx=grid.dx;  Ni=Nx-2;
    t_vec = 0:dt:params.T;  Nt=length(t_vec);

    % Grünwald-Letnikov weights: w(1)=1, w(j)=w(j-1)*(1-(1+alpha)/(j-1))
    gl = zeros(Nt,1);  gl(1)=1;
    for j=2:Nt, gl(j)=gl(j-1)*(1-(1+alpha)/(j-1)); end

    e  = ones(Ni,1);
    L  = spdiags([e,-2*e,e],[-1,0,1],Ni,Ni)/dx^2;
    mu = dt^(-alpha);

    % System: (mu*I - c*L)*u^n = -mu*history + f
    A = mu*speye(Ni) - c*L;

    u = params.ic(x(:));
    u = bc_handler(u, params.bc, grid, 0);
    U = zeros(Nx,Nt);  U(:,1)=u;
    u_hist = zeros(Ni,Nt);  u_hist(:,1)=u(2:end-1);

    for n=1:Nt-1
        mem = zeros(Ni,1);
        for j=1:n, mem=mem+gl(j+1)*u_hist(:,n+1-j); end

        fv = fsrc(params.f, grid.xi, t_vec(n+1));
        rhs = -mu*mem + fv;

        if strcmpi(params.bc.type,'dirichlet')
            uL=sv(params.bc,'left', t_vec(n+1));
            uR=sv(params.bc,'right',t_vec(n+1));
            rhs(1)  =rhs(1)  +c*uL/dx^2;
            rhs(end)=rhs(end)+c*uR/dx^2;
        end

        u(2:end-1) = A\rhs;
        u = bc_handler(u,params.bc,grid,t_vec(n+1));
        U(:,n+1)=u;
        u_hist(:,n+1)=u(2:end-1);
    end

    sol.u=U; sol.x=x; sol.t=t_vec;
    sol.info.alpha=alpha; sol.info.scheme='GL-L1 implicit';
end

function fv=fsrc(f,xv,t)
    if isequal(f,0)||(isnumeric(f)&&f==0), fv=zeros(numel(xv),1); return; end
    try fv=f(xv(:),t); catch, fv=zeros(numel(xv),1); end
end
function val=sv(bc,s,t)
    if isfield(bc,s),v=bc.(s);if isa(v,'function_handle'),val=v(t);else val=v;end
    else,val=0;end
end
