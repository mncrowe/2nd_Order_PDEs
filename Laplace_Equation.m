% Solve the Laplace Equation:
%
% u_xx + u_yy = 0
%
% u(0, y) = u(1, y) = 0
% u(x, 0) = f(x)
% u_y(x, 1) = g(x)

addpath("functions")

% Set f(x) and g(x):

f = @(x) exp(-50*(x-0.5).^2);
g = @(x) 2*sin(pi*x);

% Create numerical grid:

Nx = 64;
Ny = 64;

[Mx, x] = grid_spectral(2, Nx, [0, 1]);
[My, y] = grid_spectral(2, Ny, [0, 1]);

% Create boundary condition u(0, y) = 0:

BC1_x = create_BC(1, Nx, 1);
BC1_x_val = zeros(Ny, 1);

% Create boundary condition u(1, y) = 0:

BC2_x = create_BC(1, Nx, Nx);
BC2_x_val = zeros(Ny, 1);

% Create boundary condition u(x, 0) = f(x):

BC1_y = create_BC(1, Ny, 1);
BC1_y_val = f(x);

% Create boundary condition u_y(x, 1) = g(x):

BC2_y = create_BC(2, Ny, Ny, My);
BC2_y_val = g(x);

% Define linear left-hand side operator d^2/dx^2 + d^2/dy^2:

D = create_operator(2, eye(Nx), 0, My^2, 0) + ...
    create_operator(2, Mx^2, 0, eye(Ny), 0);

D = apply_BC(D, [Nx, Ny], [BC1_x; BC2_x], [BC1_y; BC2_y]);

% Create right-hand side term 0 (this includes BC values):

R = create_RHS(zeros(Nx, Ny), [BC1_x_val, BC2_x_val], ...
    [BC1_y_val, BC2_y_val]);

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
