R = 1000;               
C = 100e-6;            
V0 = 10;                
tfinal = 1.0;           
dt_list = [0.1, 0.05, 0.01, 0.001]; 

% Preallocation
max_errors_charge = zeros(size(dt_list));
max_errors_discharge = zeros(size(dt_list));


figure('Position', [100, 100, 1000, 800]);

%% --- SIMULATION LOOP ---
for i = 1:length(dt_list)
    dt = dt_list(i);
    
    % --- CHARGING LOOP ---
    t = 0;
    V = 0; % Starts at 0V
    t_charge = [];
    V_charge = [];
    
    while t <= tfinal
        % Need to change to dynamic allocation here to simplify
        % implementation for multiple timesteps. Originally, A matrix with
        % row size equal to the 1/minimum time step was used.
        t_charge(end+1) = t;
        V_charge(end+1) = V;
        
        dVdt_charge = (V0 - V) / (R * C);
        V = V + dt * dVdt_charge;
        t = t + dt;
    end
    
    % --- DISCHARGING LOOP ---
    t = 0;
    V = V0; % Starts fully charged at V0 = 10V
    t_discharge = [];
    V_discharge = [];
    
    while t <= tfinal
        t_discharge(end+1) = t; 
        V_discharge(end+1) = V;
        
        dVdt_discharge = -V / (R * C); % V0 is 0 during discharge
        V = V + dt * dVdt_discharge;
        t = t + dt;
    end
    
    % ANALYTICAL SOLNS
    V_a_charge = V0 * (1 - exp(-t_charge / (R * C)));
    V_a_discharge = V0 * exp(-t_discharge / (R * C));
    
    % Relative Errors (eps prevents division by zero at V=0)
    err_charge = abs(V_charge - V_a_charge) ./ (abs(V_a_charge) + eps) .* 100;
    err_discharge = abs(V_discharge - V_a_discharge) ./ (abs(V_a_discharge) + eps) .* 100;
    
    % Store maximum errors
    max_errors_charge(i) = max(err_charge);
    max_errors_discharge(i) = max(err_discharge);
    
    % NUMERICAL DATA
    % Charging
    subplot(3, 2, 1);
    plot(t_charge, V_charge, '--.', 'DisplayName', ['dt = ' num2str(dt)]);
    hold on;
    
    % Discharging
    subplot(3, 2, 2);
    plot(t_discharge, V_discharge, '--.', 'DisplayName', ['dt = ' num2str(dt)]);
    hold on;
    
    % Error vs Time - Charging
    subplot(3, 2, 3);
    semilogy(t_charge, err_charge, '-', 'DisplayName', ['Charge dt % = ' num2str(dt)]);
    hold on;

    % Error vs Time - Discharging
    subplot(3,2,4)
    semilogy(t_discharge, err_discharge, '--', 'DisplayName', ['Discharge dt %= ' num2str(dt)]);
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
loglog(dt_list, max_errors_charge, 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Max Charge Error');
hold on;
loglog(dt_list, max_errors_discharge, 'bs--', 'LineWidth', 1.5, 'DisplayName', 'Max Discharge Error');
title('Maximum Relative Error vs. Timestep');
xlabel('dt (s)');
ylabel('Max Relative Error');
grid on; 
legend('Location', 'southoutside');

