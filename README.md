# uDESolve_PDE
Universal PDE Solver — extension of the uDESolve engine. Covers parabolic, hyperbolic, elliptic, nonlinear, and fractional PDEs in 1D, 2D, and 3D
# uDESolve_PDE

Universal PDE Solver — extension of the uDESolve engine.  
Covers parabolic, hyperbolic, elliptic, nonlinear, and fractional PDEs in 1D, 2D, and 3D.

## Supported PDEs

| Type | Equation |
|---|---|
| Heat | `u_t = c·∇²u + f` |
| Wave | `u_tt = c²·∇²u` |
| Laplace | `∇²u = 0` |
| Poisson | `∇²u = f(x)` |
| Burgers | `u_t + u·u_x = ν·u_xx` |
| Fractional | `D^α u = c·∇²u + f,  α∈(0,1)` |

## Schemes

| Scheme | Details |
|---|---|
| `fdm` | FTCS, Crank-Nicolson, ADI (2D), LOD (3D), Leapfrog, Upwind, Lax-Wendroff, Fromm, Grünwald-Letnikov |
| `fem` | Linear (1D), bilinear Q1 (2D), Crank-Nicolson, Newmark-β |
| `spectral` | Chebyshev pseudospectral + RK4 / CN |
| `etdrk4` | Cox-Matthews ETDRK4 for stiff and fractional problems |
| `analytical` | Fourier series, d'Alembert, separation of variables |

## Boundary Conditions

- **Dirichlet** — direct enforcement  
- **Neumann** — 2nd-order ghost-point  
- **Robin** — `α·u + β·∂u/∂n = γ`, ghost-point elimination  

All BC values accept scalars or time-dependent function handles `@(t)`.

## Quick Start

```matlab
% 1D Heat equation — Crank-Nicolson
p.scheme     = 'fdm';
p.dim        = 1;
p.domain     = [0 1];
p.N          = 100;
p.T          = 0.5;
p.dt         = 1e-3;
p.ic         = @(x) sin(pi*x);
p.bc.type    = 'dirichlet';
p.bc.left    = 0;
p.bc.right   = 0;
p.coeff.c    = 0.01;
p.verbose    = true;

sol = pde_solver('heat', p);

% 2D Poisson — FEM
p2.scheme     = 'fem';
p2.dim        = 2;
p2.domain     = [0 1; 0 1];
p2.N          = 30;
p2.f          = @(x,y) -2*(x.*(1-x) + y.*(1-y));
p2.bc.type    = 'dirichlet';
p2.bc.left    = 0;  p2.bc.right  = 0;
p2.bc.bottom  = 0;  p2.bc.top    = 0;

sol2 = pde_solver('poisson', p2);

% 1D Burgers — Chebyshev spectral
p3.scheme    = 'spectral';
p3.dim       = 1;
p3.domain    = [-1 1];
p3.N         = 64;
p3.T         = 1;
p3.dt        = 1e-3;
p3.ic        = @(x) -sin(pi*x);
p3.bc.type   = 'dirichlet';
p3.bc.left   = 0;  p3.bc.right = 0;
p3.coeff.nu  = 0.1/pi;

sol3 = pde_solver('burgers', p3);
```

## Structure
uDESolve_PDE/
├── pde_solver.m
├── core/
│   ├── pde_validator.m      # 18-check validation layer
│   ├── pde_parser.m
│   ├── grid_generator.m
│   └── bc_handler.m
└── schemes/
    ├── fdm/
    │   ├── fdm_heat.m
    │   ├── fdm_wave.m
    │   ├── fdm_elliptic.m
    │   ├── fdm_burgers.m
    │   └── fdm_fractional.m
    ├── fem/
    │   ├── fem_assembler.m
    │   ├── fem_elliptic.m
    │   ├── fem_heat.m
    │   └── fem_wave.m
    └── spectral/
        ├── cheb_matrix.m
        ├── spectral_heat.m
        ├── spectral_elliptic.m
        └── spectral_burgers.m
```

## Author

Abdulrasheed Kamal  
University of Lagos, Nigeria.  
[github.com/AbdulrasheedKamal](https://github.com/AbdulrasheedKamal)
