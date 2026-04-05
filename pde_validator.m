function result = pde_validator(pde_type, params)
% Validates pde_type and params before any solver runs.

    VALID_TYPES   = {'heat','wave','laplace','poisson','burgers','fractional'};
    VALID_SCHEMES = {'fdm','fem','spectral','etdrk4','analytical'};
    TIME_DEP      = {'heat','wave','burgers','fractional'};
    RESERVED      = {'i','j','pi','eps','inf','nan','true','false','ans'};

    result.valid     = true;
    result.error_msg = '';
    result.warnings  = {};

    % 1. pde_type: non-empty string, known value
    if ~(ischar(pde_type) || isstring(pde_type)) || isempty(strtrim(char(pde_type)))
        result = fail(result, 'pde_type must be a non-empty string.');
        return
    end
    pde_type = lower(strtrim(char(pde_type)));
    if ~ismember(pde_type, VALID_TYPES)
        result = fail(result, sprintf('Unknown PDE type "%s". Valid: %s.', ...
            pde_type, strjoin(VALID_TYPES, ' | ')));
        return
    end

    % 2. params is a struct
    if ~isstruct(params)
        result = fail(result, 'params must be a struct.');
        return
    end

    %  scheme
    if ~isfield(params, 'scheme')
        result = fail(result, 'params.scheme required (fdm | fem | spectral | etdrk4 | analytical).');
        return
    end
    if ~ismember(lower(params.scheme), VALID_SCHEMES)
        result = fail(result, sprintf('Unknown scheme "%s".', params.scheme));
        return
    end

    % 4. if domain exists and is numeric
    if ~isfield(params, 'domain')
        result = fail(result, 'params.domain required.');
        return
    end
    if ~isnumeric(params.domain) || isempty(params.domain)
        result = fail(result, 'params.domain must be a non-empty numeric array.');
        return
    end

    % 5. domain is finite
    if any(~isfinite(params.domain(:)))
        result = fail(result, 'params.domain contains Inf or NaN.');
        return
    end

    % 6. dim is 1, 2, or 3
    dim = 1;
    if isfield(params, 'dim'), dim = params.dim; end
    if ~isnumeric(dim) || ~isscalar(dim) || ~ismember(dim, [1 2 3])
        result = fail(result, 'params.dim must be 1, 2, or 3.');
        return
    end

    % 7. domain shape matches dim
    dom = params.domain;
    if dim == 1 && numel(dom) ~= 2
        result = fail(result, '1D domain must be [xmin xmax].');
        return
    end
    if dim == 2 && ~(isequal(size(dom),[2 2]) || (isvector(dom) && numel(dom)==4))
        result = fail(result, '2D domain must be [xmin xmax; ymin ymax].');
        return
    end
    if dim == 3 && ~isequal(size(dom),[3 2])
        result = fail(result, '3D domain must be a 3x2 matrix [xmin xmax; ymin ymax; zmin zmax].');
        return
    end

    % 8. each dimension: lower < upper
    dom_r = reshape(dom, [], 2);
    for d = 1:size(dom_r,1)
        if dom_r(d,1) >= dom_r(d,2)
            result = fail(result, sprintf('Domain dim %d: %.4g >= %.4g (lower bound must be < upper).', ...
                d, dom_r(d,1), dom_r(d,2)));
            return
        end
    end

    % 9. N: integer >= 4
    if isfield(params, 'N')
        N = params.N;
        if ~isnumeric(N) || ~isfinite(N) || any(N(:) ~= floor(N(:))) || any(N(:) < 4)
            result = fail(result, 'params.N must be a finite integer >= 4.');
            return
        end
    end

    % 10. BC struct and type
    if ~isfield(params, 'bc') || ~isstruct(params.bc)
        result = fail(result, 'params.bc struct required.');
        return
    end
    if ~isfield(params.bc, 'type')
        result = fail(result, 'params.bc.type required (dirichlet | neumann | robin).');
        return
    end
    bc_type = lower(params.bc.type);
    if ~ismember(bc_type, {'dirichlet','neumann','robin'})
        result = fail(result, sprintf('Invalid BC type "%s".', params.bc.type));
        return
    end

    % 11. Contradictory condition: pure Neumann on elliptic PDE is ill-posed
    if ismember(pde_type, {'laplace','poisson'}) && strcmp(bc_type, 'neumann')
        result = fail(result, ...
            'Pure Neumann BCs on an elliptic PDE are ill-posed (solution not unique). Use Dirichlet or Robin.');
        return
    end

    % 12. Robin BCs need alpha, beta — both finite, not both zero
    if strcmp(bc_type, 'robin')
        if ~isfield(params.bc,'alpha') || ~isfield(params.bc,'beta')
            result = fail(result, 'Robin BCs require params.bc.alpha and params.bc.beta.');
            return
        end
        if ~isfinite(params.bc.alpha) || ~isfinite(params.bc.beta)
            result = fail(result, 'params.bc.alpha and params.bc.beta must be finite scalars.');
            return
        end
        if params.bc.alpha == 0 && params.bc.beta == 0
            result = fail(result, 'Robin BC: alpha=0 and beta=0 is degenerate — no condition is imposed.');
            return
        end
    end

    % 13. Equation-order vs IC count
    % Wave is 2nd order in time: ic2 (du/dt at t=0) must be a function handle if provided.
    if strcmp(pde_type, 'wave') && isfield(params, 'ic2')
        if ~isa(params.ic2, 'function_handle')
            result = fail(result, 'params.ic2 (initial velocity for wave eq.) must be a function handle.');
            return
        end
    end
    % Fractional order α must be strictly in (0,1) — not 0 (ODE), not 1 (classic PDE)
    if strcmp(pde_type, 'fractional')
        if ~isfield(params,'coeff') || ~isfield(params.coeff,'alpha_frac')
            result = fail(result, 'Fractional PDE requires params.coeff.alpha_frac in (0,1).');
            return
        end
        af = params.coeff.alpha_frac;
        if ~isscalar(af) || ~isfinite(af) || af <= 0 || af >= 1
            result = fail(result, sprintf('params.coeff.alpha_frac = %.4g must be in the open interval (0,1).', af));
            return
        end
    end

    % 14. Time-dependent PDEs: T, dt, ic
    if ismember(pde_type, TIME_DEP)
        if ~isfield(params,'T') || ~isfield(params,'dt')
            result = fail(result, 'Time-dependent PDEs require params.T and params.dt.');
            return
        end
        if ~isscalar(params.T) || ~isfinite(params.T) || params.T <= 0
            result = fail(result, 'params.T must be a positive finite scalar.');
            return
        end
        if ~isscalar(params.dt) || ~isfinite(params.dt) || params.dt <= 0
            result = fail(result, 'params.dt must be a positive finite scalar.');
            return
        end
        if params.dt >= params.T
            result = fail(result, sprintf('params.dt (%.4g) must be strictly less than params.T (%.4g).', ...
                params.dt, params.T));
            return
        end
        nsteps = ceil(params.T / params.dt);
        if nsteps > 1e6
            result.warnings{end+1} = sprintf('%.0e time steps — expect long runtime.', nsteps);
        end
        if ~isfield(params, 'ic')
            result = fail(result, 'Time-dependent PDEs require params.ic (function handle).');
            return
        end
        if ~isa(params.ic, 'function_handle')
            result = fail(result, 'params.ic must be a function handle, e.g. @(x) sin(pi*x).');
            return
        end
    end

    % 15. coeff fields: finite scalars only
    if isfield(params, 'coeff')
        co = params.coeff;
        for fn = {'c','nu','alpha_frac'}
            f = fn{1};
            if isfield(co, f)
                v = co.(f);
                if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
                    result = fail(result, sprintf('params.coeff.%s must be a finite scalar.', f));
                    return
                end
            end
        end
        if isfield(co,'c') && co.c < 0
            result.warnings{end+1} = 'params.coeff.c < 0: negative diffusivity — likely ill-posed.';
        end
        if isfield(co,'nu') && co.nu < 0
            result = fail(result, 'params.coeff.nu (viscosity) must be non-negative.');
            return
        end
    end

    % 16. source function type
    if isfield(params, 'f') && ~isequal(params.f, 0)
        if ~isa(params.f, 'function_handle') && ~(isnumeric(params.f) && isscalar(params.f))
            result = fail(result, 'params.f must be a function handle or the scalar 0.');
            return
        end
    end

    % 17. Reserved variable names in params fields
    fnames = fieldnames(params);
    for k = 1:numel(fnames)
        if ismember(lower(fnames{k}), RESERVED)
            result = fail(result, sprintf('params.%s clashes with a reserved MATLAB name.', fnames{k}));
            return
        end
    end

    % 18. Scheme-PDE compatibility
    scheme = lower(params.scheme);
    if strcmp(pde_type,'burgers') && strcmp(scheme,'fem')
        result = fail(result, 'FEM is not implemented for Burgers. Use fdm or spectral.');
        return
    end
    if strcmp(pde_type,'fractional') && ~ismember(scheme,{'fdm','etdrk4'})
        result = fail(result, sprintf('Fractional PDEs only support fdm or etdrk4, not "%s".', scheme));
        return
    end
    if strcmp(pde_type,'wave') && strcmp(scheme,'spectral')
        result.warnings{end+1} = 'Spectral wave solver not implemented; falling back is not automatic — use fdm or fem.';
    end

    if isfield(params,'verbose') && params.verbose && ~isempty(result.warnings)
        for w = 1:numel(result.warnings)
            fprintf('[validator] %s\n', result.warnings{w});
        end
    end
end

function result = fail(result, msg)
    result.valid     = false;
    result.error_msg = msg;
end
