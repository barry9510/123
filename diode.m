global Is n V_T R;  %用global來共享參數
Is = 1e-12;   
n = 1.5;      
V_T = 25.85e-3;  
R = 1000;    

function I = newton_diode(Vs)
    global Is n V_T R; 
    I = 1e-6;  % 初始猜測值
    tol = 1e-9;  % 收斂條件
    max_iter = 100;  % 最大迭代次數

    for iter = 1:max_iter
        f_I = diode_eq(I, Vs);  % 計算 f(I)
        df_I = diode_eq_derivative(I, Vs);  % 計算 f'(I)
        I_new = I - f_I / df_I;  % 牛頓法更新

        if abs(I_new - I) < tol  % 若收斂則返回結果
            I = I_new;
            return;
        end
        
        I = I_new;  % 更新 I 值
    end
end


function f = diode_eq(I, Vs)
    global Is n V_T R;
    f = I - Is * (exp((Vs - I * R) / (n * V_T)) - 1);
end


function df = diode_eq_derivative(I, Vs)
    global Is n V_T R;
    df = 1 + (Is * R / (n * V_T)) * exp((Vs - I * R) / (n * V_T)); 
end


Vs_values = linspace(0, 1, 100);  % 設定電壓範圍 (0V ~ 0.8V)
I_values = arrayfun(@newton_diode, Vs_values); %arrayfun是用來對陣列的每個元素套用指定的函數 存進I_values這個陣列中
Vd_values = Vs_values - I_values * R; % 計算二極體電壓 V_D

disp(" Vs (V)      I (mA)      V_D (V) ");
disp([Vs_values' I_values' * 1e3 Vd_values']); 

figure;
plot(Vs_values, I_values * 1e3, 'b', 'LineWidth', 2); % 轉換成 mA
xlabel('V_s (V)');
ylabel('I (mA)');
title('求解 I-V 曲線');
grid on;
legend('數值解'); 
