%global Is n V_T R L C;
Is = 1e-12;    % 反向飽和電流 (A)
n = 1.5;       % 理想因子
V_T = 25.85e-3; % 熱電壓 (V)
R = 300;       % 電阻 (Ω)
L = 1e-2;     % 電感 (H) 
C = 1e-6;      % 電容 (F)

% 時間設定
dt = 1e-6;  % 時間步長 (s)
T = 3e-3;   % 總模擬時間 (s)
N = round(T / dt);  % 四捨五入步數
time = linspace(0, T, N); % 時間軸

% 初始條件
Vs = 0.4;      % 直流電壓源
I = zeros(1, N);  % 電流 I
Vd = zeros(1, N); % 二極體電壓 V_D
Ic = zeros(1, N); 


for k = 1:N-1
    Id = Is * (exp(Vd(k) / (n * V_T)) - 1); % 二極體電流
    dI_dt = (Vs - I(k)*R - Vd(k)) / L; % 電流變化率
    dVd_dt = (I(k) - Id) / C; % 二極體電壓變化率
    Ic(k) = I(k)-Id;
    
    I(k+1) = I(k) + dI_dt * dt;
    Vd(k+1) = Vd(k) + dVd_dt * dt;
end

figure;
plot(time * 1e3, I * 1e3, 'b', 'LineWidth', 2);
xlabel('時間 (ms)');
ylabel('電流 I (mA)');
title('瞬態電流');
grid on;

figure;
plot(time * 1e3, Vd, 'r', 'LineWidth', 2);
xlabel('時間 (ms)');
ylabel('二極體電壓 V_D (V)');
title('RLC 電路中二極體的瞬態電壓');
grid on;


figure;
plot(time * 1e3, (Ic) * 1e3, 'b', 'LineWidth', 2);
xlabel('時間 (ms)');
ylabel('電流 I-Id (mA)');
title('電容電流');
grid on;


