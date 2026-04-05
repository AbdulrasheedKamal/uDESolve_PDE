function sol = fdm_elliptic(pde_info, params)

    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core'));

    dim  = params.dim;
    grid = grid_generator(params);

    if dim == 1
        x=grid.x; Nx=grid.Nx; dx=grid.dx; Ni=Nx-2;
        e=ones(Ni,1);
        L=spdiags([e,-2*e,e],[-1,0,1],Ni,Ni)/dx^2;
        fv=fsrc(params.f, grid.xi, 0);
        rhs=fv;
        if strcmpi(params.bc.type,'dirichlet')
            rhs(1)  =rhs(1)  -sv(params.bc,'left', 0)/dx^2;
            rhs(end)=rhs(end)-sv(params.bc,'right',0)/dx^2;
        end
        u=zeros(Nx,1); u(2:end-1)=L\rhs;
        u=bc_handler(u,params.bc,grid,0);
        sol.u=u; sol.x=x;

    elseif dim == 2
        x=grid.x; Nx=grid.Nx; dx=grid.dx;
        y=grid.y; Ny=grid.Ny; dy=grid.dy;
        Ni=Nx-2; Nj=Ny-2; N=Ni*Nj;
        ex=ones(Ni,1); ey=ones(Nj,1);
        Dx=spdiags([ex,-2*ex,ex],[-1,0,1],Ni,Ni)/dx^2;
        Dy=spdiags([ey,-2*ey,ey],[-1,0,1],Nj,Nj)/dy^2;
        L2=kron(speye(Nj),Dx)+kron(Dy,speye(Ni));

        [Xi,Yi]=meshgrid(grid.xi,grid.yi);
        if isequal(params.f,0)||(isnumeric(params.f)&&params.f==0)
            rhs=zeros(N,1);
        else
            rhs=reshape(params.f(Xi,Yi),N,1);
        end

        if strcmpi(params.bc.type,'dirichlet')
            uL=sv(params.bc,'left',0); uR=sv(params.bc,'right',0);
            uB=sv(params.bc,'bottom',0); uT=sv(params.bc,'top',0);
            for j=1:Nj
                b=(j-1)*Ni;
                rhs(b+1) =rhs(b+1) -uL/dx^2;
                rhs(b+Ni)=rhs(b+Ni)-uR/dx^2;
            end
            rhs(1:Ni)=rhs(1:Ni)-uB/dy^2;
            rhs((Nj-1)*Ni+1:N)=rhs((Nj-1)*Ni+1:N)-uT/dy^2;
        end

        u=zeros(Ny,Nx);
        u(2:end-1,2:end-1)=reshape(L2\rhs,Nj,Ni);
        u=bc_handler(u,params.bc,grid,0);
        sol.u=u; sol.x=x; sol.y=y;

    elseif dim == 3
        x=grid.x; Nx=grid.Nx; dx=grid.dx;
        y=grid.y; Ny=grid.Ny; dy=grid.dy;
        z=grid.z; Nz=grid.Nz; dz=grid.dz;
        Ni=Nx-2; Nj=Ny-2; Nk=Nz-2; N=Ni*Nj*Nk;
        ex=ones(Ni,1); ey=ones(Nj,1); ez=ones(Nk,1);
        Dx=spdiags([ex,-2*ex,ex],[-1,0,1],Ni,Ni)/dx^2;
        Dy=spdiags([ey,-2*ey,ey],[-1,0,1],Nj,Nj)/dy^2;
        Dz=spdiags([ez,-2*ez,ez],[-1,0,1],Nk,Nk)/dz^2;
        L3=kron(speye(Nk),kron(speye(Nj),Dx)) ...
          +kron(speye(Nk),kron(Dy,speye(Ni))) ...
          +kron(kron(Dz,speye(Nj)),speye(Ni));

        [Xi,Yi,Zi]=meshgrid(grid.xi,grid.yi,grid.zi);
        if isequal(params.f,0)||(isnumeric(params.f)&&params.f==0)
            rhs=zeros(N,1);
        else
            rhs=reshape(params.f(Xi,Yi,Zi),N,1);
        end

        if strcmpi(params.bc.type,'dirichlet')
            uL=sv(params.bc,'left',0); uR=sv(params.bc,'right',0);
            uB=sv(params.bc,'bottom',0); uT=sv(params.bc,'top',0);
            uF=sv(params.bc,'front',0); uBk=sv(params.bc,'back',0);
            for k=1:Nk, for j=1:Nj
                b=(k-1)*Nj*Ni+(j-1)*Ni;
                rhs(b+1) =rhs(b+1) -uL/dx^2;
                rhs(b+Ni)=rhs(b+Ni)-uR/dx^2;
            end, end
            for k=1:Nk, for i=1:Ni
                rhs((k-1)*Nj*Ni+i)          =rhs((k-1)*Nj*Ni+i)          -uB/dy^2;
                rhs((k-1)*Nj*Ni+(Nj-1)*Ni+i)=rhs((k-1)*Nj*Ni+(Nj-1)*Ni+i)-uT/dy^2;
            end, end
            for j=1:Nj, for i=1:Ni
                rhs((j-1)*Ni+i)          =rhs((j-1)*Ni+i)          -uF/dz^2;
                rhs((Nk-1)*Nj*Ni+(j-1)*Ni+i)=rhs((Nk-1)*Nj*Ni+(j-1)*Ni+i)-uBk/dz^2;
            end, end
        end

        u=zeros(Ny,Nx,Nz);
        u(2:end-1,2:end-1,2:end-1)=reshape(L3\rhs,Nj,Ni,Nk);
        u=bc_handler(u,params.bc,grid,0);
        sol.u=u; sol.x=x; sol.y=y; sol.z=z;
    end
end

function fv=fsrc(f,xv,t)
    if isequal(f,0)||(isnumeric(f)&&f==0), fv=zeros(numel(xv),1); return; end
    try fv=f(xv(:),t); catch, fv=zeros(numel(xv),1); end
end
function val=sv(bc,s,t)
    if isfield(bc,s),v=bc.(s);if isa(v,'function_handle'),val=v(t);else val=v;end
    else,val=0;end
end
