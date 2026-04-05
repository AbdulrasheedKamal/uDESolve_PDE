function sol = fdm_heat(pde_info, params)
% FDM solver for: u_t = c*nabla^2(u) + f
%
% sub_scheme:
%   'explicit'  FTCS,             stable if r = c*dt/dx^2 <= 0.5
%   'cn'        Crank-Nicolson,   unconditionally stable  [default]
%   'adi'       Peaceman-Rachford ADI (2D)
%   'lod'       Locally 1D / Douglas splitting (3D)

    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core'));

    dim  = params.dim;
    grid = grid_generator(params);
    c    = pde_info.coeff.c;
    dt   = params.dt;
    sub  = lower(params.sub_scheme);

    t_vec = 0:dt:params.T;
    Nt    = length(t_vec);

    if dim == 1
        x  = grid.x;  Nx = grid.Nx;  dx = grid.dx;
        r  = c*dt/dx^2;

        if strcmp(sub,'explicit') && r > 0.5
            warning('uDESolve:FDM:Heat', 'FTCS: r=%.4f > 0.5, unstable.', r);
        end

        u = params.ic(x(:));
        u = bc_handler(u, params.bc, grid, 0);
        U = zeros(Nx, Nt);
        U(:,1) = u;

        if strcmp(sub,'explicit')
            for n = 1:Nt-1
                u_new = u;
                u_new(2:end-1) = u(2:end-1) ...
                    + r*(u(3:end) - 2*u(2:end-1) + u(1:end-2)) ...
                    + dt*fsrc(params.f, x(2:end-1), t_vec(n));
                u = bc_handler(u_new, params.bc, grid, t_vec(n+1));
                U(:,n+1) = u;
            end
        else
            % Crank-Nicolson
            Ni = Nx-2;
            e  = ones(Ni,1);
            A  = spdiags([-r/2*e, (1+r)*e, -r/2*e], [-1,0,1], Ni, Ni);
            B  = spdiags([ r/2*e, (1-r)*e,  r/2*e], [-1,0,1], Ni, Ni);
            for n = 1:Nt-1
                fn  = fsrc(params.f, x(2:end-1), t_vec(n));
                fn1 = fsrc(params.f, x(2:end-1), t_vec(n+1));
                rhs = B*u(2:end-1) + (dt/2)*(fn+fn1);
                if strcmpi(params.bc.type,'dirichlet')
                    rhs(1)   = rhs(1)   + (r/2)*sv(params.bc,'left',  t_vec(n+1));
                    rhs(end) = rhs(end) + (r/2)*sv(params.bc,'right', t_vec(n+1));
                end
                u(2:end-1) = A \ rhs;
                u = bc_handler(u, params.bc, grid, t_vec(n+1));
                U(:,n+1) = u;
            end
        end

        sol.u = U;  sol.x = x;  sol.t = t_vec;
        sol.info.r = r;

    elseif dim == 2
        % Peaceman-Rachford ADI
        x=grid.x; Nx=grid.Nx; dx=grid.dx;
        y=grid.y; Ny=grid.Ny; dy=grid.dy;
        rx=c*dt/dx^2;  ry=c*dt/dy^2;
        Nxi=Nx-2; Nyi=Ny-2;
        ex=ones(Nxi,1); ey=ones(Nyi,1);

        Ax=spdiags([-rx/2*ex,(1+rx)*ex,-rx/2*ex],[-1,0,1],Nxi,Nxi);
        Bx=spdiags([ rx/2*ex,(1-rx)*ex, rx/2*ex],[-1,0,1],Nxi,Nxi);
        Ay=spdiags([-ry/2*ey,(1+ry)*ey,-ry/2*ey],[-1,0,1],Nyi,Nyi);
        By=spdiags([ ry/2*ey,(1-ry)*ey, ry/2*ey],[-1,0,1],Nyi,Nyi);
        [Lax,Uax]=lu(Ax); [Lay,Uay]=lu(Ay);

        [X,Y]=meshgrid(x,y);
        u=params.ic(X,Y);
        u=bc_handler(u,params.bc,grid,0);
        U_store=zeros(Ny,Nx,Nt); U_store(:,:,1)=u;

        for n=1:Nt-1
            th=t_vec(n)+dt/2; tn1=t_vec(n+1);
            u_half=zeros(Ny,Nx);
            for j=2:Ny-1
                row  = u(j,2:end-1)';
                expl = ry/2*(u(j-1,2:end-1)'-2*u(j,2:end-1)'+u(j+1,2:end-1)');
                rhs  = Bx*row+expl;
                rhs(1)  =rhs(1)  +(rx/2)*sv(params.bc,'left', th);
                rhs(end)=rhs(end)+(rx/2)*sv(params.bc,'right',th);
                u_half(j,2:end-1)=(Uax\(Lax\rhs))';
            end
            u_half=bc_handler(u_half,params.bc,grid,th);
            u_new=zeros(Ny,Nx);
            for i=2:Nx-1
                col  = u_half(2:end-1,i);
                expl = rx/2*(u_half(2:end-1,i-1)-2*u_half(2:end-1,i)+u_half(2:end-1,i+1));
                rhs  = By*col+expl;
                rhs(1)  =rhs(1)  +(ry/2)*sv(params.bc,'bottom',tn1);
                rhs(end)=rhs(end)+(ry/2)*sv(params.bc,'top',   tn1);
                u_new(2:end-1,i)=Uay\(Lay\rhs);
            end
            u=bc_handler(u_new,params.bc,grid,tn1);
            U_store(:,:,n+1)=u;
        end

        sol.u=U_store; sol.x=x; sol.y=y; sol.t=t_vec;
        sol.info.rx=rx; sol.info.ry=ry;

    elseif dim == 3
        % LOD (Douglas splitting)
        x=grid.x; Nx=grid.Nx; dx=grid.dx;
        y=grid.y; Ny=grid.Ny; dy=grid.dy;
        z=grid.z; Nz=grid.Nz; dz=grid.dz;
        Nxi=Nx-2; Nyi=Ny-2; Nzi=Nz-2;
        ex=ones(Nxi,1); ey=ones(Nyi,1); ez=ones(Nzi,1);
        Ax=spdiags([-c*dt/(3*dx^2)*ex,(1+2*c*dt/(3*dx^2))*ex,-c*dt/(3*dx^2)*ex],[-1,0,1],Nxi,Nxi);
        Ay=spdiags([-c*dt/(3*dy^2)*ey,(1+2*c*dt/(3*dy^2))*ey,-c*dt/(3*dy^2)*ey],[-1,0,1],Nyi,Nyi);
        Az=spdiags([-c*dt/(3*dz^2)*ez,(1+2*c*dt/(3*dz^2))*ez,-c*dt/(3*dz^2)*ez],[-1,0,1],Nzi,Nzi);
        [Lx,Ux]=lu(Ax); [Ly,Uy]=lu(Ay); [Lz,Uz]=lu(Az);

        [X,Y,Z]=meshgrid(x,y,z);
        u=params.ic(X,Y,Z);
        U_store=zeros(Ny,Nx,Nz,Nt); U_store(:,:,:,1)=u;

        for n=1:Nt-1
            u1=u;
            for k=2:Nz-1, for j=2:Ny-1
                r=squeeze(u(j,2:end-1,k))'; u1(j,2:end-1,k)=(Ux\(Lx\r))';
            end, end
            u2=u1;
            for k=2:Nz-1, for i=2:Nx-1
                col=squeeze(u1(2:end-1,i,k)); u2(2:end-1,i,k)=Uy\(Ly\col);
            end, end
            u3=u2;
            for j=2:Ny-1, for i=2:Nx-1
                tube=squeeze(u2(j,i,2:end-1)); u3(j,i,2:end-1)=Uz\(Lz\tube);
            end, end
            u=bc_handler(u3,params.bc,grid,t_vec(n+1));
            U_store(:,:,:,n+1)=u;
        end

        sol.u=U_store; sol.x=x; sol.y=y; sol.z=z; sol.t=t_vec;
    end
end

function fv = fsrc(f, xv, t)
    if isequal(f,0)||(isnumeric(f)&&f==0), fv=zeros(size(xv(:))); return; end
    try fv=f(xv(:),t); catch, fv=zeros(size(xv(:))); end
end

function val = sv(bc,s,t)
    if isfield(bc,s), v=bc.(s); if isa(v,'function_handle'),val=v(t);else val=v;end
    else, val=0; end
end
