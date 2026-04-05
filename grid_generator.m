function grid = grid_generator(params)

    dim = params.dim;
    N   = params.N;
    dom = params.domain;

    dom_r = reshape(dom, [], 2);

    Ni = N(1);
    Nj = N(min(2,end));
    Nk = N(min(3,end));

    grid.x  = linspace(dom_r(1,1), dom_r(1,2), Ni+2);
    grid.dx = grid.x(2) - grid.x(1);
    grid.Nx = Ni+2;
    grid.Ni = Ni;
    grid.xi = grid.x(2:end-1);
    grid.dim = dim;

    if dim >= 2
        grid.y  = linspace(dom_r(2,1), dom_r(2,2), Nj+2);
        grid.dy = grid.y(2) - grid.y(1);
        grid.Ny = Nj+2;
        grid.Nj = Nj;
        grid.yi = grid.y(2:end-1);
    end

    if dim == 3
        grid.z  = linspace(dom_r(3,1), dom_r(3,2), Nk+2);
        grid.dz = grid.z(2) - grid.z(1);
        grid.Nz = Nk+2;
        grid.Nk = Nk;
        grid.zi = grid.z(2:end-1);
    end

    if dim == 2
        [grid.X,  grid.Y]  = meshgrid(grid.x,  grid.y);
        [grid.Xi, grid.Yi] = meshgrid(grid.xi, grid.yi);
    elseif dim == 3
        [grid.X,  grid.Y,  grid.Z]  = meshgrid(grid.x,  grid.y,  grid.z);
        [grid.Xi, grid.Yi, grid.Zi] = meshgrid(grid.xi, grid.yi, grid.zi);
    end

    grid.h_min = grid.dx;
    if dim >= 2, grid.h_min = min(grid.h_min, grid.dy); end
    if dim == 3, grid.h_min = min(grid.h_min, grid.dz); end
end
