% CEE/MAE M20 Summer session A 2020
% Cassandra Cantu       UID: 305-100-205
%% Homework 1
% Function: let the user pick which problem in hw 1 they want to execute
% Pick 1: calculate oblate spheroid surface area from inputs
% Pick 2: identify cell's neighbors from inputs

% Clears command window
clc; close all; clear all;
%% User choice
prob_num = input('Which problem do you want to run? Enter 1 for oblate spheroid calculations or 2 for neighbor identification. \n');
switch prob_num
%% Code for 1. oblate spheroid calculations
    case 1
         r1 = input('Please enter the equatorial radius of the oblate spheroid: ');
         %checks if r1 is a non-zero & non-negative number
         if (r1<=0)
            fprintf('The surface area cannot be solved given this value. Please enter an appropiate value.')
            return
         end
         %checks if r2 is a non-zero & non-negative number
         r2 = input('Please enter the polar radius of the oblate spheroid: ');
         if (r2<=0)
            fprintf('The surface area cannot be solved given this value. Please enter an appropiate value.')
            %return
         end
         %checks if r2 is less than r1 for calculation to work
         if (r2<r1)
         %define gamma function
             g = acos(r2/r1);
         %define surface area function
             exact = 2*pi*(r1^2+((r2^2/sin(g))*log(cos(g)/(1-sin(g)))));
         %define surface area approximation function
            approx = 4*pi*((r1+r2)/2)^2;
         %prints results
            fprintf('The exact surface area is: %10.4e \n', exact);
            fprintf('The surface area approximation is: %10.4e \n', approx)
         %when condition r2<r1 is not met
         else 
             disp('The surface area cannot be solved given these values. Please enter appropriate values.')
            return
         end
%% Code for 2. neighbor identification
    case 2
        M = input('Please enter the number of rows in the array: ');
        %checks if M is greater than 2 & is an integer  
        if (M<2)|| M-floor(M)~=0
                 error('Invalid value for number of rows')
                 return
        end
        N = input('Please enter the number of columns in the array: ');
        %checks if N is greater than 2 & is an integer  
        if (N<2)||N-floor(N)~=0
                error('Invalid value for number of columns')
                return
        end
        P = input('Please enter the number of the cell ID in the array: ');
        %checks if P is in range of valid cell # & is an integer
        if (P<=1)||(P>=(M*N))|| P-floor(P)~=0
                error('Invalid value for number of target cell')
                return
        end
   %calculates all possible neighbors of P (regardless of location)
        %upper left diagonal neighbor
        P_upleft = P-M-1;
        %left neighbor
        P_left = P-M;
        %lower left diagonal neighbor
        P_downleft = P-M+1;
        %above neighbor
        P_up = P-1;
        %below neighbor
        P_down = P+1;
        %upper right diagonal neighbor
        P_upright = P+M-1;
        %right neighbor
        P_right = P+M;
        %lower right diagonal neighbor
        P_downright = P+M+1;
   %Identify location of P
    %corners (3 neighrbors total)
      %P at upper left corner
        if (P==1)
            neighbors = [P_down, P_right, P_downright];
      %P at lower left corner
        elseif (P==M)
            neighbors = [P_up, P_upright, P_right];
      %P at upper right corner
        elseif (P==M*N-(M-1))
            neighbors = [P_left, P_downleft, P_down];
      %P at lower right corner
        elseif (P==M*N)
            neighbors = [P_upleft, P_left, P_up];
    %edges (5 neighbors total)
      %P at left edge
        elseif (1<P)&&(P<M)
            neighbors = [P_up, P_down, P_upright, P_right, P_downright];
      %P at upper edge
        elseif (1<P<(M*N-(M-1)))&&(mod(P,M)==1)
            neighbors = [P_left, P_downleft, P_down, P_right, P_downright];
      %P at lower edge
        elseif (M<P<N*M)&&(mod(P,M)==0)
            neighbors = [P_upleft, P_left, P_up, P_upright, P_right];
      %P at right edge
        elseif (M*N-(M-1)<P)&&(P<M*N)
            neighbors = [P_upleft, P_left, P_downleft, P_up, P_down];
    %P is in interior (8 neighbors total)
        else
            neighbors = [P_upleft, P_left, P_downleft, P_up, P_down, P_upright, P_right, P_downright];
        end
    %Prints cell ID & neighbors in command window
        fprintf('Cell ID:   %d', P)
        fprintf('\n')
        fprintf('Neighbors: ')
        fprintf('%d  ', neighbors)
        fprintf('\n')
%checks if user doesn't enter 1 or 2 for problem number
    otherwise
        fprintf('Incorrect problem number input. Please enter a different number. \n');
end