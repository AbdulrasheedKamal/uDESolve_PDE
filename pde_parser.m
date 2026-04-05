function info = pde_parser(pde_type, params)

    info.type  = lower(pde_type);
    info.dim   = params.dim;
    info.coeff = params.coeff;

    switch info.type
        case 'heat'
            info.class       = 'parabolic';
            info.is_timedep  = true;
            info.order_time  = 1;
            info.order_space = 2;
            info.eqn         = 'u_t = c * nabla^2(u) + f';

        case 'wave'
            info.class       = 'hyperbolic';
            info.is_timedep  = true;
            info.order_time  = 2;
            info.order_space = 2;
            info.eqn         = 'u_tt = c^2 * nabla^2(u)';
            if isfield(params,'ic2')
                info.ic2 = params.ic2;
            else
                info.ic2 = @(varargin) zeros(size(varargin{1}));
            end

        case 'laplace'
            info.class       = 'elliptic';
            info.is_timedep  = false;
            info.order_time  = 0;
            info.order_space = 2;
            info.eqn         = 'nabla^2(u) = 0';

        case 'poisson'
            info.class       = 'elliptic';
            info.is_timedep  = false;
            info.order_time  = 0;
            info.order_space = 2;
            info.eqn         = 'nabla^2(u) = f(x)';

        case 'burgers'
            info.class       = 'nonlinear';
            info.is_timedep  = true;
            info.order_time  = 1;
            info.order_space = 2;
            info.eqn         = 'u_t + u*u_x = nu*u_xx';

        case 'fractional'
            info.class       = 'fractional';
            info.is_timedep  = true;
            info.order_time  = params.coeff.alpha_frac;
            info.order_space = 2;
            info.eqn         = 'D^alpha(u) = c*nabla^2(u) + f, alpha in (0,1)';

        otherwise
            error('uDESolve:Parser', 'Unknown PDE type "%s".', pde_type);
    end

    info.has_source = isfield(params,'f') && ~isequal(params.f,0) ...
                      && ~(isnumeric(params.f) && params.f == 0);

    if isfield(params,'verbose') && params.verbose
        fprintf('[parser] %s  |  %s\n', info.class, info.eqn);
    end
end
