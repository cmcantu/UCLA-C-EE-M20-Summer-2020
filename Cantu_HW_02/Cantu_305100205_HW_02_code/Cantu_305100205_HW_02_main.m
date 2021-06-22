% CEE/MAE M20 Summer session A 2020
% Cassandra Cantu       UID: 305-100-205
%% Homework 2
% Function: let the user pick which problem in hw 2 they want to execute
% Pick 1: three species predator-prey problem
% Pick 2: pocket change problem

% Clears command window
clc; close all; clear all;

%% User Choice
prob_num = input('Which problem do you want to run? Enter 1 for the three species problem or 2 for the pocket change problem: ');
switch prob_num
%% Code for 1. Three species problem
    case 1
format longG;
% set initial conditions, time conditions, & iteration count
        x_0 = 2;    y_0 = 2.49;     z_0 = 1.5;
        delta_t = 0.005;    t_final = 12;
        iter = 0;
%assign initial conditions
    x_old = x_0;
    y_old = y_0;
    z_old = z_0;
%prints heading for output and first line @t=0 (inital conditions)
    fprintf('Time\t X\t Y\t Z');
    fprintf('\n')
    fprintf('0.0\t');
    fprintf('%3.2f\t', x_old, y_old, z_old);
    fprintf('\n')
%loop over time, starting at 0, step size as 0.005, end at t=12)
tic
    for t = 0:delta_t:t_final
        %forward euler method
        x_new = x_old + ((0.75*x_old*(1-(x_old/20))-(1.5*x_old*y_old)-(0.5*x_old*z_old))*delta_t);
        y_new = y_old + (y_old*(1-(y_old/25))-(0.75*x_old*y_old)-(1.25*y_old*z_old))*delta_t;
        z_new = z_old + (1.5*z_old*(1-(z_old/30))-(x_old*z_old)-(y_old*z_old))*delta_t;
%update variables
        x_old = x_new;
        y_old = y_new;
        z_old = z_new;
        iter = iter + 1;
%prints time-varying populations for delta_t = 0.5 s (every 100 iterations)
        if (mod(iter, 100)==0)
                fprintf('%2.1f\t', t);
                fprintf('%3.2f\t', x_new);
                fprintf('%3.2f\t', y_new); 
                fprintf('%3.2f\t', z_new);
                fprintf('\n');
        end
    end
    Time = toc;
%% Code for 2. Pocket change problem
    case 2
%Initialize total sum of coins
        sum = 0;
%loop to calculate # of coins for each starting amount of change (1-99 cents)
        for start_change = 0:1:99
         change = start_change;
%initalize variables, restarts for each iteration
         Q = 0; D = 0; N = 0; P = 0;
%checks if rest of change is > 0
            while change > 0
%counts number of quarters until rest of change<25
                if change - 25 >= 0
                    Q = Q+1;
                    change = change - 25;
%counts number of dimes until rest of change<10
                 elseif change - 10 >= 0
                     D = D+1;
                     change = change - 10;
% counts number of nickels until rest of change<5
                 elseif change - 5 >= 0
                     N = N+1;
                     change = change - 4;
%counts number of pennies until rest of change<1
                elseif change - 1 >= 0
                    P=P+1;
                    change = change - 1;
                end
%calculates average number of coins you expect to receive in change
            end
             total_coins = Q + D + N + P;
             sum = sum + total_coins;
            coin_mean = sum/100;
        end
           fprintf('Average Number of Coins = %.2f \n', coin_mean);
%checks if user doesn't enter 1 or 2 for problem number
    otherwise
        fprintf('Incorrect problem number input. Please enter a different number. \n');
end