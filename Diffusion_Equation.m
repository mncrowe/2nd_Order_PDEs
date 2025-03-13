% Solve the Diffusion Equation:
%
% u_t - a * u_xx = 0
%
% u(0, y) = u(1, y) = 0
% u(x, 0) = f(x)

addpath("functions")

% Set diffusivity and f(x):

a = 0.05;
f = @(x) exp(-50*(x-0.5).^2);

% Create numerical grid:

Nx = 64;
Ny = 32;

[Mx, x] = grid_spectral(2, Nx, [0, 1]);
[My, y] = grid_spectral(2, Ny, [0, 1]);

% Create boundary condition u(0, y) = 0:

BC1_x = create_BC(1, Nx, 1);
BC1_x_val = zeros(Ny, 1);

% Create boundary condition u(1, y) = 0:

BC2_x = create_BC(1, Nx, Nx);
BC2_x_val = zeros(Ny, 1);

% Create initial condition u(x, 0) = f(x):

BC_y_val = f(x);

% Define linear left-hand side operator d/dy - a * d^2/dx^2:

D = create_operator(2, eye(Nx), 0, My, 0) - ...
    a * create_operator(2, Mx^2, 0, eye(Ny), 0);

D = apply_BC(D, [Nx, Ny], [BC1_x; BC2_x], 0);
D(1:Nx, :) = [eye(Nx), zeros(Nx, Nx*(Ny-1))];

% Create right-hand side term 0 (this includes BC values):

R = create_RHS(zeros(Nx, Ny), [BC1_x_val, BC2_x_val], 0);
R(1:Nx) = BC_y_val;

% Solve equation:

u = D \ R;
u = reshape(u, Nx, Ny);

% Plot solution:

close all
surf(x, y, u');
xlabel('x'); ylabel('y'); zlabel('u');
view(0, 90); axis equal
xlim([x(1) x(end)]); ylim([y(1) y(end)])
colormap(cmap2()); colorbar
