% CEE/MAE M20 Summer session A 2020
% Cassandra Cantu       UID: 305-100-205
%% Homework 6
% Function: let the user pick which problem in hw 6 they want to execute
% Pick 1: game of life
% Pick 2: Euler-Bernoulli beam bending

% Clears command window
clc; close all; clear all;
%% User choice
prob_num = input('Which problem do you want to run? Enter:\n1 for the game of life or\n2 for Euler-Bernoulli beam bending \n');
switch prob_num
%% Code for 1. game of life
    case 1
   %set number of rows
        x_grid = 150;
   %set number of columns
        y_grid = 200;
   %initialize lifeGrid (current) and newGrid (nex gen)
        lifeGrid = rand(x_grid, y_grid); %0 to 1
        newGrid = zeros(x_grid, y_grid);
   %go through cells, assign 0 or 1 depending on probability
        for i=1:1:x_grid
            for j=1:1:y_grid
              %alive cells
                if lifeGrid(i,j)>0.9
                    lifeGrid(i,j) = 1;
              %dead cells
                else
                    lifeGrid(i,j) = 0;
                end
            end
        end
        
   %initial condition plot
        figure(1)
        imagesc(lifeGrid)
        title('Game of Life Initial Conditions')
        
   %define number of generations
        Max_iter = 300;
   %initialize iteration count
        iter = 0;
   %initialize empty array to store number of live cells per iter
       live_cells = zeros(Max_iter, 1);
        while iter<Max_iter
          %update iteration
            iter = iter+1;
         %go through cells
           %rows
            for i=1:1:x_grid
               %check boundary conditions
                if i==1 %1st row
                    im1=x_grid;
                else
                    im1=i-1;
                end
                
                if i==x_grid %last row
                    ip1=1;
                else
                    ip1=i+1;
                end
           %columns
            for j=1:1:y_grid
               %check boundary conditions
                if j==1 %1st column
                    jm1=y_grid;
                else
                    jm1= j-1;
                end
                
                if j==y_grid %last row
                    jp1=1;
                else
                    jp1=j+1;
                end
             
        %calculate living neighbors (%upper left, upper, and upper right
                                     %left and right
                                    %lower left, lower, and lower right)
            live_ngb = lifeGrid(im1,jm1)+lifeGrid(im1,j)+lifeGrid(im1,jp1)+...
                    lifeGrid(i, jm1)+lifeGrid(i, jp1)+...
                    lifeGrid(ip1, jm1)+lifeGrid(ip1, j)+lifeGrid(ip1,jp1);
                    
         %check live or dead
             if lifeGrid(i,j)==1 %originally alive
               %live cell lives when 2 or 3 living neighbors
                 if live_ngb==2|| live_ngb==3
                     newGrid(i,j)=1;
               %live cell dies when <2 or >3 living neighbors
                 else
                     newGrid(i,j)=0;
                 end
              else %originally dead
                %dead cell becomes alive when 3 living neighbors
                  if live_ngb==3
                      newGrid(i,j)=1;
                  end
             end
            end
            end
            
          %display animation
            figure(2)
            imagesc(newGrid)
            drawnow
            title(['Game of Life (generation = ' num2str(iter) ')'], 'FontSize', 24)
            xlabel('Number of Columns', 'FontSize', 20)
            ylabel('Number of Rows', 'FontSize', 20)
            
          %update grid after 1 generation
            lifeGrid = newGrid;
          %store number of live cells
            live_cells(iter) = sum(newGrid, "all");
        end
          %plot number of living cells as function of time
            figure(3)
            p = plot(live_cells, 'LineWidth', 2);
            xlabel('Time (generation)', 'FontSize', 20)
            ylabel('Number of Living Cells', 'FontSize', 20)
            title('Number of Living Cells over Time', 'FontSize', 24)
            
                    
%% Code for 2. Euler-Bernoulli beam bending
    case 2
    %initialize dimensions
        R = 0.013; %outer radius
        r = 0.011; %inner radius
        L = 1; %beam length
        E = 70e9; %Young's modulus
        d = 0.75; %where applied force is placed
        p = -2000; %applied force
        %calculate moment of inertia
        I = (pi/4)*(R^4-r^4);
        
    %set number of nodes
        n_grid = 20;
        %initialize node coordinate
        x = zeros(n_grid,1);
        dx = L/(n_grid-1);
        for i = 2:1:n_grid
            x(i) = x(i-1)+dx;
        end
        
     %create A matrix
        A = zeros(n_grid, n_grid);
       %boundary conditions
        A(1,1) = 1;
        A(n_grid,n_grid) = 1;
       %band, A as coefficients
        for i = 2:1:n_grid-1
            A(i,i-1) = 1;
            A(i,i) = -2;
            A(i,i+1) = 1;
        end
        
        %define right-hand side
        M = zeros(n_grid,1);
        for i=1:1:n_grid
            if x(i) <= d
                M(i) = -p*(L-d)*x(i)/L;
            else
                M(i) = -p*d*(L-x(i))/L;
            end
        end
        
        %Create b vector
        b = dx^2*M/(E*I);
        
        %solve y (deflection)
        y = A\b;
        
        %plot
        figure(1)
        title('Euler-Bernoulli beam bending')
      %deflection
        subplot(2,1,1) 
        title('Deflection')
        plot(x,y,'o-','LineWidth',3)
        xlabel('x [m]')
        ylabel('Deflection [m]')
      %moment
        subplot(2,1,2);
        title('Moment')
        plot(x,M,'b--', 'LineWidth', 3)
        xlabel('x [m]')
        ylabel('Moment [Nm]')
    %calculate error
      %maximum displacement using min function
        max_ydis = min(y);
      %initialize location of max displacement
        loc=0;
      %find & store location
        for i=1:1:n_grid
            if y(i) == max_ydis
                loc = x(i);
            end
        end
     %print results
        fprintf('Calculated maximum displacement: %f m\n', max_ydis);
        fprintf('Location of maximum displacement from left side: %f m\n', loc);
      %calculate maximum displacement w/ theoretical model
        c = min(d, L-d);
        theo_ymax = p*c*((L^2-c^2)^1.5)/(9*sqrt(3)*E*I*L);
        fprintf('Theoretical maximum displacement: %f m\n', theo_ymax)
       %calculate percent error 
        per_err = abs((max_ydis-theo_ymax)/theo_ymax)*100;
        fprintf('Percent error of calculation method compared to theoretical solution: ');
        fprintf('%f', per_err);
        fprintf('%%\n')
    otherwise
        fprintf('Incorrect problem number input. Please enter a different number. \n');
end