% CEE/MAE M20 Summer session A 2020
% Cassandra Cantu       UID: 305-100-205
%% Homework 5
% Function: let the user pick which problem in hw 5 they want to execute
% Pick 1: shared birthday problem
% Pick 2: random walk collisions

% Clears command window
clc; close all; clear all;
%% User choice
prob_num = input('Which problem do you want to run? Enter:\n1 for the shared birthday problem or\n2 for random walk collision \n');
switch prob_num
%% Code for 1. shared birthday problem
    case 1
    n_trials=1e4; %number of trials
    n_people = zeros(n_trials, 1); %empty array to store # people for each trial
    for trial = 1:1:n_trials
      %initialization
        match_flag = 0;
        iter = 1; %number of people in group
        iterMax = 1e4;
        group = ceil(rand*365);
       %loop for generating new birthdays
        while match_flag == 0 && iter<iterMax
          %add person to group
            iter = iter+1;
          %generate new birthday
            nb = ceil(rand*365);
          %search all birthdays
                for k = 1:1:length(group)
                  %look for 2 in same week
                    if abs(nb - group(k))<7 || abs(nb - group(k) - 365)<7
                        match_flag = 1;
                    end
                end
           %update the birthday
             new_group = [group; nb];
             group = new_group;
        end
      %save total number of people for current trial
        n_people(trial)=iter;
    end
  %calculate median
    med_val = ceil(median(n_people));
    fprintf('Median Number of People = %02d\n', med_val);
  %print histogram
        figure(1)
        histogram(n_people)
        grid on
        title('Possibility of Shared Birthdays within a Week', 'Fontsize', 24)
        xlabel('Number of People in Group', 'Fontsize', 20)
        ylabel('Number of Trials', 'Fontsize', 20)
        set(gcf,'Position',[100 100 1500 750])
        set(gca, 'LineWidth', 3, 'FontSize', 20)
%% Code for 2. random walk collisions
    case 2
        %initialization
        n_trials = 5000;
        n_iter = zeros(n_trials, 1);
        BC = [5, -5, -5, 5]; %[up, down, left, right]
        %loop through trials to find collision
        for k=1:1:n_trials
            %initialize A&B walkers x & y positions
            Ax = -5; Ax_i = Ax;
            Ay = 0; Ay_i = Ay;
            Bx = 5; Bx_i = Bx;
            By = 0; By_i = By;
            %clear figure for iteration
            %%note: comment out bc not displaying window for 5000 iter
           % clf;
          %set max number of steps for walkers
            n_steps = 1000;
          %initialize condition marker & step count
            collision_flag = 0;
            i=0;
            %use function RandWalk_2D to get new position
            while collision_flag == 0 && i<n_steps
                [Ax_ip1, Ay_ip1] = RandWalk_2D(Ax_i,Ay_i,BC);
                [Bx_ip1, By_ip1] = RandWalk_2D(Bx_i,By_i,BC);
                
           %create a grid for A (step i)
            Ax_ival = [Ax_i - 0.5, Ax_i + 0.5, Ax_i + 0.5, Ax_i - 0.5];
            Ay_ival = [Ay_i - 0.5, Ay_i - 0.5, Ay_i + 0.5, Ay_i + 0.5];
           %create a grid for A (step i+1)
            Ax_ip1val = [Ax_ip1 - 0.5, Ax_ip1 + 0.5, Ax_ip1 + 0.5, Ax_ip1 - 0.5];
            Ay_ip1val = [Ay_ip1 - 0.5, Ay_ip1 - 0.5, Ay_ip1 + 0.5, Ay_ip1 + 0.5];
            
           %create a grid for B (step i)
            Bx_ival = [Bx_i - 0.5, Bx_i + 0.5, Bx_i + 0.5, Bx_i - 0.5];
            By_ival = [By_i - 0.5, By_i - 0.5, By_i + 0.5, By_i + 0.5];
           %create a grid for B (step i+1)
            Bx_ip1val = [Bx_ip1 - 0.5, Bx_ip1 + 0.5, Bx_ip1 + 0.5, Bx_ip1 - 0.5];
            By_ip1val = [By_ip1 - 0.5, By_ip1 - 0.5, By_ip1 + 0.5, By_ip1 + 0.5];
           
           %% create plot 
           %%note: comment section out so running 5000 trials will not take so long
%             figure(1)
%             hold on
%             set(gca, 'xtick', -5:1:5)
%             set(gca, 'ytick', -5:1:5)
%             grid on
%             xlim([-5.5, 5.5])
%             ylim([-5.5, 5.5])
%             axis square
% %            %fill in current step of A
%             fill(Ax_ival, Ay_ival, 'r')
% %            %fill in i+1 step of A
%             fill(Ax_ip1val, Ay_ip1val, 'b')
% %            %fill in current step of B
%             fill(Bx_ival, By_ival, 'y')
% %            %fill in i+1 step of B
%            fill(Bx_ip1val, By_ip1val, 'g')
%             title('2D Random Walk','FontSize',24)
%             set (gcf, 'Position', [100 100 1500 750])
%             set(gca, 'LineWidth', 3, 'FontSize', 20)
%             hold off
            
         %update position & iteration
            Ax_i = Ax_ip1;  Ay_i = Ay_ip1;
           Bx_i = Bx_ip1;  By_i = By_ip1;
            i = i+1;
            %check for collision & store
                if Ax_i == Bx_i && Ay_i == By_i
                    collision_flag = 1;
                    n_iter(k) = i;
                end
            end
        end
        med_val = median(n_iter);
        fprintf('Median = %2i\n', med_val);
    otherwise
          fprintf('Incorrect problem number input. Please enter a different number. \n');
end