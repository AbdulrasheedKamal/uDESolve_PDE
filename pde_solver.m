function sol = pde_solver(pde_type, params)
% Main dispatcher for uDESolve PDE extension.
%
% Usage:
%   sol = pde_solver('heat', params)
%
% PDE types : heat | wave | laplace | poisson | burgers | fractional
% Schemes   : fdm  | fem  | spectral | etdrk4  | analytical
%
% Required params fields:
%   .scheme, .domain, .dim, .bc  (always)
%   .N, .T, .dt, .ic             (time-dep problems)
%   .coeff.c, .coeff.nu, ...     (PDE coefficients)
%
% sol fields:
%   .u, .x [.y .z], .t, .scheme, .pde_type, .info

    vcheck = pde_validator(pde_type, params);
    if ~vcheck.valid
        error('uDESolve:PDE:InvalidInput', '%s', vcheck.error_msg);
    end

    params   = set_defaults(params);
    t_start  = tic;
    pde_info = pde_parser(pde_type, params);
    scheme   = lower(params.scheme);
    ptype    = lower(pde_type);

    root = fileparts(mfilename('fullpath'));

    switch scheme
        case 'fdm'
            addpath(fullfile(root,'schemes','fdm'));
            sol = dispatch_fdm(ptype, pde_info, params);
        case 'fem'
            addpath(fullfile(root,'schemes','fem'));
            sol = dispatch_fem(ptype, pde_info, params);
        case 'spectral'
            addpath(fullfile(root,'schemes','spectral'));
            sol = dispatch_spectral(ptype, pde_info, params);
        case 'etdrk4'
            addpath(fullfile(root,'schemes','etdrk4'));
            sol = etdrk4_solver(ptype, pde_info, params);
        case 'analytical'
            addpath(fullfile(root,'schemes','analytical'));
            sol = dispatch_analytical(ptype, pde_info, params);
        otherwise
            error('uDESolve:PDE:UnknownScheme', 'Unknown scheme: %s', scheme);
    end

    sol.pde_type     = ptype;
    sol.scheme       = scheme;
    sol.info.runtime = toc(t_start);
    sol.info.dim     = params.dim;
    sol.info.N       = params.N;

    if params.verbose
        fprintf('[uDESolve] %s/%s  dim=%d  N=%d  t=%.3fs\n', ...
            upper(ptype), upper(scheme), params.dim, params.N, sol.info.runtime);
    end
end

function sol = dispatch_fdm(ptype, info, p)
    switch ptype
        case 'heat',                  sol = fdm_heat(info, p);
        case 'wave',                  sol = fdm_wave(info, p);
        case {'laplace','poisson'},   sol = fdm_elliptic(info, p);
        case 'burgers',               sol = fdm_burgers(info, p);
        case 'fractional',            sol = fdm_fractional(info, p);
        otherwise, error('uDESolve:FDM', 'FDM: no solver for "%s"', ptype);
    end
end

function sol = dispatch_fem(ptype, info, p)
    switch ptype
        case 'heat',                  sol = fem_heat(info, p);
        case 'wave',                  sol = fem_wave(info, p);
        case {'laplace','poisson'},   sol = fem_elliptic(info, p);
        otherwise, error('uDESolve:FEM', 'FEM: no solver for "%s"', ptype);
    end
end

function sol = dispatch_spectral(ptype, info, p)
    switch ptype
        case 'heat',                  sol = spectral_heat(info, p);
        case {'laplace','poisson'},   sol = spectral_elliptic(info, p);
        case 'burgers',               sol = spectral_burgers(info, p);
        otherwise, error('uDESolve:Spectral', 'Spectral: no solver for "%s"', ptype);
    end
end

function sol = dispatch_analytical(ptype, info, p)
    switch ptype
        case 'heat',                  sol = analytical_heat(info, p);
        case 'wave',                  sol = analytical_wave(info, p);
        case {'laplace','poisson'},   sol = analytical_elliptic(info, p);
        otherwise, error('uDESolve:Analytical', 'No analytical solver for "%s"', ptype);
    end
end

function p = set_defaults(p)
    if ~isfield(p,'verbose'),     p.verbose     = false; end
    if ~isfield(p,'dim'),         p.dim         = 1;     end
    if ~isfield(p,'N'),           p.N           = 100;   end
    if ~isfield(p,'dt'),          p.dt          = 1e-3;  end
    if ~isfield(p,'T'),           p.T           = 1.0;   end
    if ~isfield(p,'N_modes'),     p.N_modes     = 50;    end
    if ~isfield(p,'sub_scheme'),  p.sub_scheme  = 'cn';  end
    if ~isfield(p,'f'),           p.f           = 0;     end
    if ~isfield(p,'coeff'),       p.coeff       = struct(); end
    if ~isfield(p.coeff,'c'),          p.coeff.c          = 1.0; end
    if ~isfield(p.coeff,'nu'),         p.coeff.nu         = 0.01; end
    if ~isfield(p.coeff,'alpha_frac'), p.coeff.alpha_frac = 0.8;  end
    if ~isfield(p,'ic2'), p.ic2 = @(varargin) zeros(size(varargin{1})); end
end
