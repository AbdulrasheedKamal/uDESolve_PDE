function u = bc_handler(u, bc, grid, t)
% Applies boundary conditions to solution array u.
% Supports Dirichlet, Neumann (2nd-order ghost point), Robin.

    if nargin < 4, t = 0; end

    dim   = grid.dim;
    btype = lower(bc.type);

    gL = sv(bc,'left',  t);
    gR = sv(bc,'right', t);
    dx = grid.dx;

    if dim >= 2
        gB = sv(bc,'bottom',t);
        gT = sv(bc,'top',   t);
        dy = grid.dy;
    end
    if dim == 3
        gF  = sv(bc,'front',t);
        gBk = sv(bc,'back', t);
        dz  = grid.dz;
    end

    switch btype

        case 'dirichlet'
            if dim == 1
                u(1) = gL;  u(end) = gR;
            elseif dim == 2
                u(:,1) = gL;  u(:,end) = gR;
                u(1,:) = gB;  u(end,:) = gT;
            elseif dim == 3
                u(:,:,1) = gL;  u(:,:,end) = gR;
                u(:,1,:) = gB;  u(:,end,:) = gT;
                u(1,:,:) = gF;  u(end,:,:) = gBk;
            end

        case 'neumann'
            % ghost-point: u_ghost = u_interior ± 2h*g  (outward sign)
            if dim == 1
                u(1)   = u(3)     - 2*dx*gL;
                u(end) = u(end-2) + 2*dx*gR;
            elseif dim == 2
                u(:,1)   = u(:,3)     - 2*dx*gL;
                u(:,end) = u(:,end-2) + 2*dx*gR;
                u(1,:)   = u(3,:)     - 2*dy*gB;
                u(end,:) = u(end-2,:) + 2*dy*gT;
            elseif dim == 3
                u(:,:,1)   = u(:,:,3)     - 2*dx*gL;
                u(:,:,end) = u(:,:,end-2) + 2*dx*gR;
                u(:,1,:)   = u(:,3,:)     - 2*dy*gB;
                u(:,end,:) = u(:,end-2,:) + 2*dy*gT;
                u(1,:,:)   = u(3,:,:)     - 2*dz*gF;
                u(end,:,:) = u(end-2,:,:) + 2*dz*gBk;
            end

        case 'robin'
            % alpha*u + beta*(du/dn) = gamma
            % ghost-point elimination → u_bnd = (gamma*h + beta*u_inner) / (alpha*h + beta)
            al = bc.alpha;  bt = bc.beta;
            gm_L  = gv(bc,'gamma_left',  'gamma', t);
            gm_R  = gv(bc,'gamma_right', 'gamma', t);

            if dim == 1
                u(1)   = (gm_L*dx + bt*u(2))     / (al*dx + bt + eps);
                u(end) = (gm_R*dx + bt*u(end-1)) / (al*dx + bt + eps);
            elseif dim == 2
                gm_B = gv(bc,'gamma_bottom','gamma',t);
                gm_T = gv(bc,'gamma_top',   'gamma',t);
                u(:,1)   = (gm_L*dx + bt*u(:,2))     / (al*dx + bt + eps);
                u(:,end) = (gm_R*dx + bt*u(:,end-1)) / (al*dx + bt + eps);
                u(1,:)   = (gm_B*dy + bt*u(2,:))     / (al*dy + bt + eps);
                u(end,:) = (gm_T*dy + bt*u(end-1,:)) / (al*dy + bt + eps);
            elseif dim == 3
                gm_B  = gv(bc,'gamma_bottom','gamma',t);
                gm_T  = gv(bc,'gamma_top',   'gamma',t);
                gm_F  = gv(bc,'gamma_front', 'gamma',t);
                gm_Bk = gv(bc,'gamma_back',  'gamma',t);
                u(:,:,1)   = (gm_L*dx  + bt*u(:,:,2))     / (al*dx  + bt + eps);
                u(:,:,end) = (gm_R*dx  + bt*u(:,:,end-1)) / (al*dx  + bt + eps);
                u(:,1,:)   = (gm_B*dy  + bt*u(:,2,:))     / (al*dy  + bt + eps);
                u(:,end,:) = (gm_T*dy  + bt*u(:,end-1,:)) / (al*dy  + bt + eps);
                u(1,:,:)   = (gm_F*dz  + bt*u(2,:,:))     / (al*dz  + bt + eps);
                u(end,:,:) = (gm_Bk*dz + bt*u(end-1,:,:)) / (al*dz  + bt + eps);
            end

        otherwise
            error('uDESolve:BC', 'Unknown BC type: "%s"', btype);
    end
end

function val = sv(bc, side, t)
    if isfield(bc,side), v=bc.(side);
        if isa(v,'function_handle'), val=v(t); else val=v; end
    else, val=0; end
end

function val = gv(bc, primary, fallback, t)
    if isfield(bc,primary), val=sv(bc,primary,t);
    elseif isfield(bc,fallback), val=sv(bc,fallback,t);
    else, val=0; end
end
