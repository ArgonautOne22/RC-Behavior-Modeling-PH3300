R = 1000;               
C = 100e-6;            
V0 = 10;                
tfinal = 1.0;           
dt_list = [0.1, 0.05, 0.01, 0.001]; 

% Preallocation
max_errors_charge_fe = zeros(size(dt_list));
max_errors_discharge_fe = zeros(size(dt_list));
max_errors_charge_rk = zeros(size(dt_list));
max_errors_discharge_rk = zeros(size(dt_list));


figure('Position', [100, 100, 1000, 800]);

%% --- SIMULATION LOOP ---
for i = 1:length(dt_list)
    dt = dt_list(i);
    
    % --- CHARGING LOOP ---
    t = 0;
    V_fe = 0; % Starts at 0V
    V_rk = 0;
    t_charge = [];
    V_charge_fe = [];
    V_charge_rk = [];
    
    while t <= tfinal
        % Need to change to dynamic allocation here to simplify
        % implementation for multiple timesteps. Originally, A matrix with
        % row size equal to the 1/minimum time step was used.
        t_charge(end+1) = t;
        V_charge_fe(end+1) = V_fe;
        V_charge_rk(end+1) = V_rk;
        
        dVdt_charge = (V0 - V_fe) / (R * C);
        V_fe = V_fe + dt * dVdt_charge;
        t = t + dt;

        % RK4 propagation
        f = @(V) (V0 - V_rk) / (R * C);
        k1 = f(V_rk);
        k2 = f(V_rk + dt * k1 / 2);
        k3 = f(V_rk + dt * k2 / 2);
        k4 = f(V_rk + dt * k3);
        V_rk = V_rk + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
    
    % --- DISCHARGING LOOP ---
    t = 0;
    V_fe = V0; % Starts fully charged at V0 = 10V
    t_discharge = [];
    V_discharge_fe = [];
    V_discharge_rk = [];
    
    while t <= tfinal
        t_discharge(end+1) = t; 
        V_discharge_fe(end+1) = V_fe;
        V_discharge_rk(end+1) = V_rk;
        
        dVdt_discharge = -V_fe / (R * C); % V0 is 0 during discharge
        V_fe = V_fe + dt * dVdt_discharge;

        f = @(V) -V / (R * C);
        k1 = f(V_rk);
        k2 = f(V_rk + dt * k1 / 2);
        k3 = f(V_rk + dt * k2 / 2);
        k4 = f(V_rk + dt * k3);
        V_rk = V_rk + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;

        t = t + dt;
    end
    
    % ANALYTICAL SOLNS
    V_a_charge = V0 * (1 - exp(-t_charge / (R * C)));
    V_a_discharge = V0 * exp(-t_discharge / (R * C));
    
    % Relative Errors (eps prevents division by zero at V=0)
    err_charge_fe = abs(V_charge_fe - V_a_charge) ./ (abs(V_a_charge) + eps) .* 100;
    err_discharge_fe = abs(V_discharge_fe - V_a_discharge) ./ (abs(V_a_discharge) + eps) .* 100;
    err_charge_rk = abs(V_charge_rk - V_a_charge) ./ (abs(V_a_charge) + eps) .* 100;
    err_discharge_rk = abs(V_discharge_rk - V_a_discharge) ./ (abs(V_a_discharge) + eps) .* 100;
    
    % Store maximum errors
    max_errors_charge_fe(i) = max(err_charge_fe);
    max_errors_discharge_fe(i) = max(err_discharge_fe);
    max_errors_charge_rk(i) = max(err_charge_rk);
    max_errors_discharge_rk(i) = max(err_discharge_rk);
    
    % NUMERICAL DATA
    % Charging
    subplot(3, 2, 1);
    plot(t_charge, V_charge_fe, '--.', 'DisplayName', ['dt = ' num2str(dt)]);
    if (dt == 0.001)
        plot(t_charge, V_charge_rk, '--.', 'DisplayName', ['RK4 dt = ' num2str(dt)]);
    end
    hold on;
    
    % Discharging
    subplot(3, 2, 2);
    plot(t_discharge, V_discharge_fe, '--.', 'DisplayName', ['dt = ' num2str(dt)]);
    if (dt == 0.001)
        plot(t_discharge, V_discharge_rk, '--.', 'DisplayName', ['RK4 dt = ' num2str(dt)]);
    end
    hold on;
    
    % Error vs Time - Charging
    subplot(3, 2, 3);
    semilogy(t_charge, err_charge_fe, '-', 'DisplayName', ['Charge dt % = ' num2str(dt)]);
    if (dt == 0.001)
        semilogy(t_charge, err_charge_rk, '-', 'DisplayName', ['RK4 Charge dt % = ' num2str(dt)]);
    end
    hold on;

    % Error vs Time - Discharging
    subplot(3,2,4)
    semilogy(t_discharge, err_discharge_fe, '--', 'DisplayName', ['Discharge dt %= ' num2str(dt)]);
    if (dt == 0.001)
        semilogy(t_discharge, err_discharge_rk, '--', 'DisplayName', ['RK4 Discharge dt %= ' num2str(dt)]);
    end
    hold on;
end


t_dense = linspace(0, tfinal, 1000);
V_exact_charge = V0 * (1 - exp(-t_dense / (R * C)));
V_exact_discharge = V0 * exp(-t_dense / (R * C));


% Charging Curve
subplot(3, 2, 1);
plot(t_dense, V_exact_charge, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
title('RC Circuit Charging Curve');
xlabel('Time (s)'); 
ylabel('Voltage (V)');
grid on; 
legend('Location', 'southoutside');
ylim([0 12])

% Discharging Curve
subplot(3, 2, 2);
plot(t_dense, V_exact_discharge, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
title('RC Circuit Discharging Curve');
xlabel('Time (s)'); 
ylabel('Voltage (V)');
grid on; 
legend('Location', 'southoutside');
ylim([-1 12])

% Error vs Time Curve
subplot(3,2,3)
title('Relative Error Charging vs. Time');
xlabel('Time (s)'); ylabel('Relative Error');
grid on; 
legend('Location', 'southoutside');
hold off;

subplot(3,2,4)
title('Relative Error Discharging vs. Time');
xlabel('Time (s)'); ylabel('Relative Error');
grid on; 
legend('Location', 'southoutside');
hold off;

% Maximum Error vs dt
subplot(3,2,5)
loglog(dt_list, max_errors_charge_fe, 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Max Charge Error');
hold on;
loglog(dt_list, max_errors_charge_rk, 'rs-', 'LineWidth', 1.5, 'DisplayName', 'Max Charge Error RK4');
loglog(dt_list, max_errors_discharge_fe, 'bo--', 'LineWidth', 1.5, 'DisplayName', 'Max Discharge Error');
loglog(dt_list, max_errors_discharge_rk, 'bs--', 'LineWidth', 1.5, 'DisplayName', 'Max Discharge Error RK4');

title('Maximum Relative Error vs. Timestep');
xlabel('dt (s)');
ylabel('Max Relative Error');
grid on; 
legend('Location', 'southoutside');



