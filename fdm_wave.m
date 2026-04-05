function sol = fdm_wave(pde_info, params)
% FDM leapfrog solver for: u_tt = c^2 * nabla^2(u)
% CFL: r = c*dt/dx <= 1 (1D), rx+ry <= 1 (2D), rx+ry+rz <= 1 (3D)
% ic = u(x,0),  ic2 = du/dt(x,0)

    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core'));

    dim  = params.dim;
    grid = grid_generator(params);
    c    = pde_info.coeff.c;
    dt   = params.dt;

    t_vec = 0:dt:params.T;
    Nt    = length(t_vec);

    if dim == 1
        x=grid.x; Nx=grid.Nx; dx=grid.dx;
        r = c*dt/dx;
        if r > 1, warning('uDESolve:FDM:Wave','CFL r=%.4f > 1, unstable.', r); end

        u0 = params.ic(x(:));
        u0 = bc_handler(u0, params.bc, grid, 0);
        g0 = pde_info.ic2(x(:));

        u1 = u0;
        u1(2:end-1) = u0(2:end-1) + dt*g0(2:end-1) ...
            + (r^2/2)*(u0(3:end)-2*u0(2:end-1)+u0(1:end-2));
        u1 = bc_handler(u1, params.bc, grid, t_vec(min(2,Nt)));

        U = zeros(Nx,Nt);  U(:,1) = u0;
        if Nt >= 2, U(:,2) = u1; end
        u_prev = u0;  u_curr = u1;

        for n = 2:Nt-1
            u_next = zeros(Nx,1);
            u_next(2:end-1) = 2*u_curr(2:end-1) - u_prev(2:end-1) ...
                + r^2*(u_curr(3:end)-2*u_curr(2:end-1)+u_curr(1:end-2));
            u_next = bc_handler(u_next, params.bc, grid, t_vec(n+1));
            u_prev = u_curr;  u_curr = u_next;
            U(:,n+1) = u_curr;
        end

        sol.u=U; sol.x=x; sol.t=t_vec;
        sol.info.CFL = r;

    elseif dim == 2
        x=grid.x; Nx=grid.Nx; dx=grid.dx;
        y=grid.y; Ny=grid.Ny; dy=grid.dy;
        rx=(c*dt/dx)^2; ry=(c*dt/dy)^2;
        if rx+ry > 1, warning('uDESolve:FDM:Wave2D','2D CFL: rx+ry=%.4f > 1.', rx+ry); end

        [X,Y]=meshgrid(x,y);
        u0=params.ic(X,Y); u0=bc_handler(u0,params.bc,grid,0);
        g0=pde_info.ic2(X,Y);
        u1=u0;
        u1(2:end-1,2:end-1)=u0(2:end-1,2:end-1)+dt*g0(2:end-1,2:end-1) ...
            +(rx/2)*(u0(2:end-1,3:end)-2*u0(2:end-1,2:end-1)+u0(2:end-1,1:end-2)) ...
            +(ry/2)*(u0(3:end,2:end-1)-2*u0(2:end-1,2:end-1)+u0(1:end-2,2:end-1));
        u1=bc_handler(u1,params.bc,grid,t_vec(min(2,Nt)));

        U_store=zeros(Ny,Nx,Nt); U_store(:,:,1)=u0;
        if Nt>=2, U_store(:,:,2)=u1; end
        u_prev=u0; u_curr=u1;

        for n=2:Nt-1
            u_next=2*u_curr-u_prev;
            u_next(2:end-1,2:end-1)=u_next(2:end-1,2:end-1) ...
                +rx*(u_curr(2:end-1,3:end)-2*u_curr(2:end-1,2:end-1)+u_curr(2:end-1,1:end-2)) ...
                +ry*(u_curr(3:end,2:end-1)-2*u_curr(2:end-1,2:end-1)+u_curr(1:end-2,2:end-1));
            u_next=bc_handler(u_next,params.bc,grid,t_vec(n+1));
            u_prev=u_curr; u_curr=u_next;
            U_store(:,:,n+1)=u_curr;
        end

        sol.u=U_store; sol.x=x; sol.y=y; sol.t=t_vec;
        sol.info.CFL=rx+ry;

    elseif dim == 3
        x=grid.x; Nx=grid.Nx; dx=grid.dx;
        y=grid.y; Ny=grid.Ny; dy=grid.dy;
        z=grid.z; Nz=grid.Nz; dz=grid.dz;
        rx=(c*dt/dx)^2; ry=(c*dt/dy)^2; rz=(c*dt/dz)^2;
        if rx+ry+rz>1, warning('uDESolve:FDM:Wave3D','3D CFL=%.4f > 1.',rx+ry+rz); end

        [X,Y,Z]=meshgrid(x,y,z);
        u0=params.ic(X,Y,Z); g0=pde_info.ic2(X,Y,Z);
        u1=u0;
        u1(2:end-1,2:end-1,2:end-1)=u0(2:end-1,2:end-1,2:end-1)+dt*g0(2:end-1,2:end-1,2:end-1) ...
            +(rx/2)*(u0(2:end-1,3:end,2:end-1)-2*u0(2:end-1,2:end-1,2:end-1)+u0(2:end-1,1:end-2,2:end-1)) ...
            +(ry/2)*(u0(3:end,2:end-1,2:end-1)-2*u0(2:end-1,2:end-1,2:end-1)+u0(1:end-2,2:end-1,2:end-1)) ...
            +(rz/2)*(u0(2:end-1,2:end-1,3:end)-2*u0(2:end-1,2:end-1,2:end-1)+u0(2:end-1,2:end-1,1:end-2));
        u1=bc_handler(u1,params.bc,grid,t_vec(min(2,Nt)));

        U_store=zeros(Ny,Nx,Nz,Nt); U_store(:,:,:,1)=u0;
        if Nt>=2, U_store(:,:,:,2)=u1; end
        u_prev=u0; u_curr=u1;

        for n=2:Nt-1
            u_next=2*u_curr-u_prev;
            u_next(2:end-1,2:end-1,2:end-1)=u_next(2:end-1,2:end-1,2:end-1) ...
                +rx*(u_curr(2:end-1,3:end,2:end-1)-2*u_curr(2:end-1,2:end-1,2:end-1)+u_curr(2:end-1,1:end-2,2:end-1)) ...
                +ry*(u_curr(3:end,2:end-1,2:end-1)-2*u_curr(2:end-1,2:end-1,2:end-1)+u_curr(1:end-2,2:end-1,2:end-1)) ...
                +rz*(u_curr(2:end-1,2:end-1,3:end)-2*u_curr(2:end-1,2:end-1,2:end-1)+u_curr(2:end-1,2:end-1,1:end-2));
            u_next=bc_handler(u_next,params.bc,grid,t_vec(n+1));
            u_prev=u_curr; u_curr=u_next;
            U_store(:,:,:,n+1)=u_curr;
        end

        sol.u=U_store; sol.x=x; sol.y=y; sol.z=z; sol.t=t_vec;
        sol.info.CFL=rx+ry+rz;
    end
end
