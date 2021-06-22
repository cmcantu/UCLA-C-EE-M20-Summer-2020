% CEE/MAE M20 Summer session A 2020
% Cassandra Cantu    
%% Homework 3
% Function: let the user pick which problem in hw 3 they want to execute
% Pick 1: pendulum physics problem
% Pick 2: DNA analysis

% Clears command window
clc; close all; clear all;
%% User choice
prob_num = input('Which problem do you want to run? Enter 1 for pendulum physics problem or 2 for DNA analysis. \n');
switch prob_num
%% Code for 1. Pendulum physics problem
    case 1
%set initial conditions
    L = 1;    g = 9.81;   theta_0 = pi/3;
%set time conditions
    t0 = 0;  dt = 0.005;    tf = 20;   
    nt = (tf - t0)/dt;
    t = linspace(t0,tf,nt); 
%initialize arrays 
    theta = zeros(1, nt);
    theta(1) = theta_0;
    omega = zeros(1, nt);
    alpha = zeros(1, nt);
    
%Forward Euler Method
    for k = 1:1:nt-1
       %calculate angular position
        theta(k+1) = theta(k)+dt*omega(k);
       %calculate angular velocity
        omega(k+1) = omega(k) - dt*g/L*sin(theta(k));
       %calculate angular acceleration
        alpha(k+1) = -g/L*sin(theta(k+1));
    end
   %calculate total energy
    h = L*(1-cos(theta));
    tot_E = g*h + 0.5*(L*omega).^2;
    
   %plot 1. kinematic vectors vs. time (Forward Euler Method) on same axis
    figure(1)
    p = plot(t,theta, 'DisplayName', 'Position');hold on;
    v = plot(t,omega, 'DisplayName', 'Velocity');hold on;
    a = plot(t,alpha, 'DisplayName', 'Acceleration');hold off;
    title('Angular Position, Velocity, & Acceleration vs. Time (Forward Euler Method)', 'FontSize', 24)
    xlabel('Time (s)', 'FontSize', 20)
    ylabel('Magnitude','FontSize', 20)
    xlim([0 20])
    ylim([-10 10])
    set(p, 'LineWidth', 4)
    set(v, 'LineWidth', 4)
    set(a, 'LineWidth', 4)
    set (gcf, 'Position', [100 100 1500 750])
    set (gca, 'LineWidth', 3, 'FontSize', 20)
    legend
   %plot 2. Total Energy vs. Time (Forward Euler Method)
    figure (2)
     e = plot(t, tot_E, 'm');
     title('Total Energy vs. Time (Forward Euler Method)')
     xlabel('Time (s)')
     ylabel('Energy (J)')
     xlim([0 20])
     ylim([4 11])
     set(e, 'LineWidth', 4)
     set (gcf, 'Position', [100 100 1500 750])
     set (gca, 'LineWidth', 3, 'FontSize', 20)
    
%Semi-Implicit Euler Method
        for k = 1:1:nt-1
           %calculate angular velocity
            omega(k+1) = omega(k) - dt*g/L*sin(theta(k));
           %calculate angular position
            theta(k+1) = theta(k)+dt*omega(k+1); 
           %calculate angular acceleration
            alpha(k+1) = -g/L*sin(theta(k+1));
        end
   %calculate total energy
    h = L*(1-cos(theta));
    tot_E = g*h + 0.5*(L*omega).^2;
    
   %plot 3. kinematic vectors vs. time (Semi-Implicit Method) on same axis
    figure(3)
    p = plot(t,theta, 'm', 'DisplayName', 'Position');hold on;
    v = plot(t,omega, 'g', 'DisplayName', 'Velocity');hold on;
    a = plot(t,alpha, 'b', 'DisplayName', 'Acceleration');hold off;
    title('Angular Position, Velocity, & Acceleration vs. Time (Semi-Implicit Euler Method)', 'FontSize', 24)
    xlabel('Time (s)', 'FontSize', 20)
    ylabel('Magnitude','FontSize', 20)
    xlim([0 20])
    ylim([-9 9])
    set(p, 'LineWidth', 4)
    set(v, 'LineWidth', 4)
    set(a, 'LineWidth', 4)
    set (gcf, 'Position', [100 100 1500 750])
    set (gca, 'LineWidth', 3, 'FontSize', 20)
    legend
   %plot 4. Total Energy vs. Time (Semi-Implicit Euler Method)
    figure (4)
     e = plot(t, tot_E, 'b');
     title('Total Energy vs. Time (Semi-Implicit Euler Method)')
     xlabel('Time (s)')
     ylabel('Energy (J)')
     xlim([0 20])
     ylim([4.86 4.95])
     set(e, 'LineWidth', 4)
     set (gcf, 'Position', [100 100 1500 750])
     set (gca, 'LineWidth', 3, 'FontSize', 20)
%% Code for 2. DNA analysis
    case 2
%determine number of bases from dna file
        file = load('chr1_sect.mat');
        dna = file.dna;
        numBases = length(dna);
%initialize startPoint, endPoint, & number of protein segments to zero
    startPoint = 0;
    endPoint = 0;
    np = 0;
    
%set stop codon counts to 0
    TAA_count = 0;
    TAG_count = 0;
    TGA_count = 0;

%look through dna sequence
    for k = 1:3:numBases-2
%Look for start
        if startPoint == 0
            if dna(k) == 1 && dna(k+1) == 4 && dna(k+2) == 3
                startPoint = k;
               %count number of protein segments
                np = np+1;
               %store values into startpoint array
                savedStart(np) = startPoint;
            end
%Look for end after start is established
        else
            if (dna(k) == 4 && dna(k+1) == 1 && dna(k+2) == 1)
                 endPoint = k;
                %store values into endpoint array;
                 savedEnd(np) = endPoint;
                %update counter for TAA
                 TAA_count = TAA_count + 1;
                %reset startPoint to zero
                 startPoint = 0;
            elseif (dna(k) == 4 && dna(k+1) == 1 && dna(k+2) == 3)
                 endPoint = k;
                %store values into endpoint array;
                 savedEnd(np) = endPoint;
                %update counter for TAA
                  TAG_count = TAG_count + 1;
                %reset startPoint to zero
                 startPoint = 0;
            elseif dna(k) == 4 && dna(k+1) == 3 && dna(k+2) == 1
                 endPoint = k;
                %store values into endpoint array;
                 savedEnd(np) = endPoint;
                %update counter for TGA
                 TGA_count = TGA_count + 1;
                %reset startPoint to zero
                 startPoint = 0;
            end
        end
    end
%calculate lengths of protein-coding segments and stores in new array
    array = savedEnd - savedStart + 3;
%calculate total protein-coding segments
    total = length(array);
%calculate average length
    avg = mean(array);
%calculate maximum length
    min = min(array);
%calculate minimum length
    max = max(array);
%calculate percentage of DNA directly used in protein-coding process
    percent = sum(array)/numBases*100;
%prints outputs
    fprintf('Total Protein-Coding Segments: %d \n', total)
    fprintf('Average Length: %.2f bases\n', avg)
    fprintf('Maximum Length: %d bases\n', max)
    fprintf('Minimum Length: %d bases\n', min)
    fprintf('Percentage of DNA Used: %.2f%%\n', percent)
%calculate most frequently used stop codon
   if  (TAA_count > TAG_count) && (TAA_count > TGA_count)
        fprintf('The most frequently used stop codon is: TAA\n');
   elseif (TAG_count > TAA_count) && (TAG_count > TGA_count)
       fprintf('The most frequently used stop codon is: TAG\n');
   else
       fprintf('The most frequently used stop codon is: TGA\n');
   end
%calculate least frequently used stop codon
   if  (TAA_count < TAG_count) && (TAA_count < TGA_count)
        fprintf('The least frequently used stop codon is: TAA\n');
   elseif (TAG_count < TAA_count) && (TAG_count < TGA_count)
       fprintf('The least frequently used stop codon is: TAG\n');
   else
       fprintf('The least frequently used stop codon is: TGA\n');
   end
%checks if user doesn't enter 1 or 2 for problem number
   otherwise
        fprintf('Incorrect problem number input. Please enter a different number. \n');
    end
