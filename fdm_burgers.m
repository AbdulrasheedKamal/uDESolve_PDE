function sol = fdm_burgers(pde_info, params)
% FDM solver for: u_t + u*u_x = nu*u_xx  (1D viscous Burgers)
%
% sub_scheme: 'upwind' | 'lax' (Lax-Wendroff) [default] | 'fromm'

    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core'));

    grid = grid_generator(params);
    nu   = pde_info.coeff.nu;
    dt   = params.dt;
    x    = grid.x;  Nx=grid.Nx;  dx=grid.dx;
    t_vec = 0:dt:params.T;  Nt=length(t_vec);

    sub = lower(params.sub_scheme);
    if ~ismember(sub,{'upwind','lax','fromm'}), sub='lax'; end

    d = nu*dt/dx^2;
    if d > 0.5, warning('uDESolve:FDM:Burgers','Diffusion number d=%.4f > 0.5.', d); end

    u = params.ic(x(:));
    u = bc_handler(u, params.bc, grid, 0);
    U = zeros(Nx,Nt);  U(:,1)=u;

    for n=1:Nt-1
        ui = u(2:end-1);
        switch sub
            case 'upwind'
                up = max(ui,0); un = min(ui,0);
                conv = up.*(ui-u(1:end-2))/dx + un.*(u(3:end)-ui)/dx;
                diff = nu*(u(3:end)-2*ui+u(1:end-2))/dx^2;
                u(2:end-1) = ui + dt*(-conv+diff);

            case 'lax'
                f_i  = ui.^2/2;
                f_ip = u(3:end).^2/2;
                f_im = u(1:end-2).^2/2;
                conv = (f_ip-f_im)/(2*dx) ...
                     - (dt/(2*dx^2))*((ui+u(3:end)).*(f_ip-f_i) ...
                                    -(u(1:end-2)+ui).*(f_i-f_im));
                diff = nu*(u(3:end)-2*ui+u(1:end-2))/dx^2;
                u(2:end-1) = ui + dt*(-conv+diff);

            case 'fromm'
                r = zeros(size(ui));
                for k=1:numel(ui)
                    denom = u(k+2)-u(k+1);
                    if abs(denom)>1e-10, r(k)=(u(k+1)-u(k))/denom; end
                end
                phi = (r+abs(r))./(1+abs(r));
                Fiph = ui.^2/2 + (dt/(2*dx))*(1-(dt/dx)*ui).*phi.*(u(3:end).^2/2-ui.^2/2)/2;
                us = u(1:end-2);
                r2 = zeros(size(ui));
                for k=1:numel(ui)
                    denom=u(k+1)-u(k);
                    if abs(denom)>1e-10, r2(k)=(u(k)-us(k))/denom; end
                end
                phi2=(r2+abs(r2))./(1+abs(r2));
                Fimh=us.^2/2+(dt/(2*dx))*(1-(dt/dx)*us).*phi2.*(ui.^2/2-us.^2/2)/2;
                conv=(Fiph-Fimh)/dx;
                diff=nu*(u(3:end)-2*ui+u(1:end-2))/dx^2;
                u(2:end-1)=ui+dt*(-conv+diff);
        end
        u = bc_handler(u, params.bc, grid, t_vec(n+1));
        U(:,n+1) = u;
    end

    sol.u=U; sol.x=x; sol.t=t_vec;
    sol.info.d=d; sol.info.sub=sub;
end
