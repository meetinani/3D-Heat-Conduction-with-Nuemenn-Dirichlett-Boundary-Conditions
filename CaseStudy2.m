% Parameters
Lx = 0.1; % Length of the device in x-direction (m)
Ly = 0.05; % Length of the device in y-direction (m)
Lz = 0.02; % Length of the device in z-direction (m)
Nx = 30; % Number of grid points in x-direction
Ny = 20; % Number of grid points in y-direction
Nz = 10; % Number of grid points in z-direction
T_initial = 25; % Initial temperature (°C)
T_source = 100; % Temperature of heat source (°C)
alpha = 1e-4; % Thermal diffusivity (m^2/s)
dt = 0.01; % Time step (s)
t_final = 10; % Final time (s)

% Grid spacing
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
dz = Lz / (Nz - 1);

% Initialize temperature matrix
T = T_initial * ones(Ny, Nx, Nz);

% Define heat source location
source_x = round(Nx / 2);
source_y = round(Ny / 2);
source_z = round(Nz / 2);

% Main time-stepping loop for unsteady-state
num_steps = t_final / dt;

for k = 1:num_steps
    % Dirichlet boundary conditions for x-direction
    T(:, 1, :) = 150; % Left boundary
    T(:, end, :) = 150; % Right boundary
    
    % Neumann boundary conditions for y-direction
    T(1, :, :) = T(2, :, :); % Top boundary
    T(end, :, :) = T(end-1, :, :); % Bottom boundary
    
    % Dirichlet boundary conditions for z-direction
    T(:, :, 1) = 75; % Front boundary
    T(:, :, end) = 75; % Back boundary
    
    % Apply heat source
    T(source_y, source_x, source_z) = T_source;
    
    % Compute temperature at next time step using finite difference
    T_new = T;
    for i = 2:Nx-1
        for j = 2:Ny-1
            for m = 2:Nz-1
                T_new(j, i, m) = T(j, i, m) + alpha * dt * ((T(j, i+1, m) - 2*T(j, i, m) + T(j, i-1, m)) / dx^2 + ...
                                                              (T(j+1, i, m) - 2*T(j, i, m) + T(j-1, i, m)) / dy^2 + ...
                                                              (T(j, i, m+1) - 2*T(j, i, m) + T(j, i, m-1)) / dz^2);
            end
        end
    end
    
    % Update temperature matrix for next time step
    T = T_new;
end

% Plot temperature distribution along x-direction
x_slice = round(Nx / 2); % Middle slice view in x-direction
figure;
imagesc(squeeze(T(:, x_slice, :)));
colorbar;
title('Temperature Distribution along X-Direction');
xlabel('y');
ylabel('z');
axis equal;

% Plot temperature distribution along y-direction
y_slice = round(Ny / 2); % Middle slice view in y-direction
figure;
imagesc(squeeze(T(y_slice, :, :)));
colorbar;
title('Temperature Distribution along Y-Direction');
xlabel('x');
ylabel('z');
axis equal;

% Plot temperature distribution along z-direction
z_slice = round(Nz / 2); % Middle slice view in z-direction
figure;
imagesc(squeeze(T(:, :, z_slice)));
colorbar;
title('Temperature Distribution along Z-Direction');
xlabel('x');
ylabel('y');
axis equal;
