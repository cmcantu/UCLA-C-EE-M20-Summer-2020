
% CEE/MAE M20 Summer session A 2020
% Cassandra Cantu       UID: 305-100-205
%% Final Project

% Clears command window
clear all;
close all;
clc;

%initialization for trials 1-8
m = [3 4 5 8 10 3 3 3];
k = [200 50 125 25 100 47 27 100];
c = [2 45 50 35 10 40 18 3];
%number of trials
n_exp = 8;
    
%time parameters
t0 = 0;
tf = 10;
dt = 1/300;
t = t0:dt:tf;
nt = length(t);

%Spring-Damping coefficents
%natural frequency
w_n = sqrt(k./m);
%damping ratio
xi = c./(2.*m.*w_n);

%initial condition
x0 = 1;

%% Free Vibration

%Free Response - 1st Order Runge-Kutta (forward Euler)
x_RK1_free = zeros(n_exp,nt);
v_RK1_free = zeros(n_exp,nt);

%run each experiment
for i=1:1:n_exp
  %set initial condition
    x_RK1_free(i,1) = x0;
    
 %go through timesteps
    for j=1:1:nt-1
      %set free vibration (0 array since we will call on elements in func.)
        f = [0 0 0];
       
       %call function to get displacement & velocity at next time
       %step
        temp = VibrationPosition([x_RK1_free(i,j), v_RK1_free(i,j)],...
        m(i), k(i), c(i), f, dt,1);
    
        x_RK1_free(i, j+1) = temp(1);
        v_RK1_free(i, j+1) = temp(2);
    end
end

%Free Response - 2nd Order Runge-Kutta
x_RK2_free = zeros(n_exp,nt);
v_RK2_free = zeros(n_exp,nt);

%run each experiment
for i=1:1:n_exp
  %set initial condition
    x_RK2_free(i,1) = x0;
    
 %go through timesteps
    for j=1:1:nt-1
      %set free vibration
        f = [0 0 0];
       
       %call function to get displacement & velocity at next time
       %step
        temp = VibrationPosition([x_RK2_free(i,j), v_RK2_free(i,j)],...
        m(i), k(i), c(i), f, dt,2);
    
        x_RK2_free(i, j+1) = temp(1);
        v_RK2_free(i, j+1) = temp(2);
        
    end
end

%Free Response - 4th Order Runge-Kutta
x_RK4_free = zeros(n_exp,nt);
v_RK4_free = zeros(n_exp,nt);

%run each experiment
for i=1:1:n_exp
  %set initial condition
    x_RK4_free(i,1) = x0;
    
 %go through timesteps
    for j=1:1:nt-1
      %set free vibration
        f = [0 0 0];
       
       %call function to get displacement & velocity at next time
       %step
        temp = VibrationPosition([x_RK4_free(i,j), v_RK4_free(i,j)],...
        m(i), k(i), c(i), f, dt, 4);
    
        x_RK4_free(i, j+1) = temp(1);
        v_RK4_free(i, j+1) = temp(2);
    end
end

%% Force Vibration
%define forcing function
a0 = 2.5;
F = a0*sin(t/(2*pi));

%Force Response - 1st Order Runge-Kutta (Forward Euler)
x_RK1_force = zeros(n_exp, nt);
v_RK1_force = zeros(n_exp,nt);

for i=1:1:n_exp
   %set initial condition
    x_RK1_force(i,1) = x0;
    
    for j=1:1:nt-1
       %define forcing function for RK1 -> F(t)/m
        f=[a0*sin(t(j)/(2*pi))/m(i) 0 0];
       
       %call function 
        temp = VibrationPosition([x_RK1_force(i,j), v_RK1_force(i,j)],...
            m(i),k(i),c(i),f,dt,1);
       
        %store values in force vibration vector
        x_RK1_force(i,j+1) = temp(1);
        v_RK1_force(i,j+1) = temp(2);
    end
end

%Force Response - 2nd Order Runge-Kutta
x_RK2_force = zeros(n_exp, nt);
v_RK2_force = zeros(n_exp,nt);

for i=1:1:n_exp
   %set initial condition
    x_RK2_force(i,1) = x0;
    
    for j=1:1:nt-1
       %define forcing function for RK2 -> f(t), f(t+0.5dt)
         f=[a0*sin(t(j)/(2*pi))/m(i) a0*sin((t(j)+0.5*dt)/(2*pi))/m(i) 0];
       
       %call function 
        temp = VibrationPosition([x_RK2_force(i,j), v_RK2_force(i,j)],...
            m(i),k(i),c(i),f,dt,2);
       
        %store values in force vibration vector
        x_RK2_force(i,j+1) = temp(1);
        v_RK2_force(i,j+1) = temp(2);
    end
end

%Force Response - 4th Order Runge-Kutta
x_RK4_force = zeros(n_exp, nt);
v_RK4_force = zeros(n_exp,nt);

for i=1:1:n_exp
   %set initial condition
    x_RK4_force(i,1) = x0;
    
    for j=1:1:nt-1
       %define forcing function array for RK4 -> f(t), f(t+0.5dt),& f(t+dt)
        f=[a0*sin(t(j)/(2*pi))/m(i) a0*sin((t(j)+0.5*dt)/(2*pi))/m(i)...
            a0*sin((t(j)+dt)/(2*pi))/m(i)];
       
       %call function 
        temp = VibrationPosition([x_RK4_force(i,j), v_RK4_force(i,j)],...
            m(i),k(i),c(i),f,dt,4);
       
        %store values in force vibration vector
        x_RK4_force(i,j+1) = temp(1);
        v_RK4_force(i,j+1) = temp(2);
    end
end

%% Plotting
for i=1:1:n_exp
    figure(i)
    hold on
    set(gcf,'Position',[15 50 1350 775])
  %Plot homogenous response v time
    subplot(1,2,1)
    hold on
    grid on
        plot(t,x_RK1_free(i,:),'ro','DisplayName','RK-1')
        plot(t,x_RK2_free(i,:),'bo','DisplayName','RK-2')
        plot(t,x_RK4_free(i,:),'go', 'DisplayName','RK-4')
        set(gca,'LineWidth',3,'FontSize',18)
        xlim([0 10])
        xlabel('Time (s)')
        ylabel('Position (m)')
        title(['Homogeneous Response (Trial ' num2str(i) ')'], 'FontSize',22)
        legend
   %Plot inhomogenous response v time
    subplot(1,2,2)
    hold on
    grid on
        plot(t,x_RK1_force(i,:),'ro','DisplayName','RK-1')
        plot(t,x_RK2_force(i,:),'bo','DisplayName','RK-2')
        plot(t,x_RK4_force(i,:),'go','DisplayName','RK-4')
        set(gca,'LineWidth',3,'FontSize',18)
        xlim([0 10])
        xlabel('Time (s)')
        ylabel('Position (m)')
        title(['Inhomogeneous Response (Trial ' num2str(i) ')'],'FontSize',22)
       legend
end
   
%% Frequency Response
%experiment with xi using given m/k/c
%xi = c./(2.*m.*w_n); %experiment by changing this eqn

%normalized frequency
lambda = logspace(-1,1,500);
n_lambda = length(lambda);
%empty matrices for Gain and Phase
Gain = zeros(n_exp,n_lambda);
PhiG = zeros(n_exp,n_lambda);

%9 so won't overwrite other figures for trials 1-8
figure(9)
set(gcf,'Position',[15 50 1350 775])
for run = 1:1:n_exp

    for j=1:1:n_lambda
       %absolute gain
        Gain(run,j) = 1/sqrt((1-lambda(j)^2)^2 + (2*xi(run)*lambda(j))^2);
       %phase shift
        PhiG(run,j) = -atan(2*xi(run)*lambda(j)/(1-lambda(j)^2));
    end
  %plot 
    subplot(2,1,1)
        semilogx(lambda,Gain(run,:),'LineWidth',3)
        hold on
        grid on
        set(gca,'LineWidth',3,'FontSize',20)
        xlabel('Normalized Frequency')
        ylabel('Amplification Factor')
        title('Bode Plot: Function Gain')
   %plot
    subplot(2,1,2)
        semilogx(lambda,PhiG(run,:),'LineWidth',3)
        hold on
        grid on
        set(gca,'LineWidth',3,'FontSize',20)
        xlabel('Normalized Frequency')
        ylabel('Phase (Degrees)')
        title('Bode Plot: Function Phase Shift')
    legendInfo{run} = ['Run' num2str(run)];
end
legend(legendInfo)


%% Animation
%Underdamped video
%create video file
vidfile1 = VideoWriter('Underdamped.mp4','MPEG-4');
vidfile1.FrameRate = 30;
open(vidfile1);

%initialize frame
nv = ceil(length(t)/10);
%set up empty arrays for t_out, x pos (homogenous), x pos (inhomogenous)
t_out = zeros(nv,1);
x_hom = zeros(nv,1);
x_inh = zeros(nv,1);
count = 0; %initialize count
s1 = 0.25; %side length for box

figure(10)
for i=1:1:length(t)
    %plot a frame every 10 sec
    if mod(i-1,10)==0
      %update count
        count=count+1;
        t_out(count) = t(i);
        x_hom(count) = x_RK4_free(1,i); %first case is underdamped
        x_inh(count) = x_RK4_force(1,i); %underdamped
       %solve y position
        y1 = [x_hom(count)+s1, x_hom(count)-s1, x_hom(count)-s1, x_hom(count)+s1]; %clockwise
        y2 = [x_inh(count)+s1, x_inh(count)-s1, x_inh(count)-s1, x_inh(count)+s1]; %clockwise
       %square
        x = [-s1 -s1 s1 s1];
       
       %homogeneous
        subplot(1,2,1)
            fill(x,y1,'b')
            xlim([-1-s1, 1+s1])
            ylim([-1-s1, 1+s1])
            xlabel('X position [m]')
            ylabel('Y position [m]')
            title('Homogeneous Free Response')
            axis square
       %inhomogeneous
        subplot(1,2,2)
            fill(x,y2,'b')
            xlim([-1-s1, 1+s1])
            ylim([-1-s1, 1+s1])
            xlabel('X position [m]')
            ylabel('Y position [m]')
            title('Inhomogeneous Force Response')
            axis square
            
            %save frame into video
            writeVideo(vidfile1,getframe(gcf));
    end
end

close(vidfile1)

%Overdamped video
%create video file
vidfile2 = VideoWriter('Overdamped.mp4','MPEG-4');
vidfile2.FrameRate = 30;
open(vidfile2);

%initialize frame
%how many frames needed
nv = ceil(length(t)/10);
t_out = zeros(nv,1);
x_hom = zeros(nv,1);
x_inh = zeros(nv,1);
count = 0;
s1 = 0.25; %side length for box

figure(11)
for i=1:1:length(t)
    %plot a frame every 10 sec
    if mod(i-1,10)==0
      %update count
        count=count+1;
        t_out(count) = t(i);
        x_hom(count) = x_RK4_free(2,i); %second case is overdamped
        x_inh(count) = x_RK4_force(2,i); %overdamped
     
        y1 = [x_hom(count)+s1, x_hom(count)-s1, x_hom(count)-s1, x_hom(count)+s1]; %clockwise
        y2 = [x_inh(count)+s1, x_inh(count)-s1, x_inh(count)-s1, x_inh(count)+s1]; %clockwise
       %square
        x = [-s1 -s1 s1 s1];
        
        subplot(1,2,1)
            fill(x,y1,'b')
            xlim([-1-s1, 1+s1])
            ylim([-1-s1, 1+s1])
            xlabel('X position [m]')
            ylabel('Y position [m]')
            title('Homogeneous Free Response')
            axis square
       
        subplot(1,2,2)
            fill(x,y2,'b')
            xlim([-1-s1, 1+s1])
            ylim([-1-s1, 1+s1])
            xlabel('X position [m]')
            ylabel('Y position [m]')
            title('Inhomogeneous Force Response')
            axis square
            
            %save frame into video
            writeVideo(vidfile2,getframe(gcf));
    end
end

close(vidfile2)

%Critically damped video
%create video file
vidfile3 = VideoWriter('CriticallyDamped.mp4','MPEG-4');
vidfile3.FrameRate = 30;
open(vidfile3);

%initialize frame
%how many frames needed
nv = ceil(length(t)/10);
t_out = zeros(nv,1);
x_hom = zeros(nv,1);
x_inh = zeros(nv,1);
count = 0;
s1 = 0.25; %side length for box

figure(12)
for i=1:1:length(t)
    %plot a frame every 10 sec
    if mod(i-1,10)==0
      %update count and variables
        count=count+1;
        t_out(count) = t(i);
        x_hom(count) = x_RK4_free(3,i); %third case is critically damped
        x_inh(count) = x_RK4_force(3,i); %critically damped
        
        y1 = [x_hom(count)+s1, x_hom(count)-s1, x_hom(count)-s1, x_hom(count)+s1]; %clockwise
        y2 = [x_inh(count)+s1, x_inh(count)-s1, x_inh(count)-s1, x_inh(count)+s1]; %clockwise
       %square
        x = [-s1 -s1 s1 s1];
        
        subplot(1,2,1)
            fill(x,y1,'b')
            xlim([-1-s1, 1+s1])
            ylim([-1-s1, 1+s1])
            xlabel('X position [m]')
            ylabel('Y position [m]')
            title('Homogeneous Free Response')
            axis square
       
        subplot(1,2,2)
            fill(x,y2,'b')
            xlim([-1-s1, 1+s1])
            ylim([-1-s1, 1+s1])
            xlabel('X position [m]')
            ylabel('Y position [m]')
            title('Inhomogeneous Force Response')
            axis square
            
            %save frame into video
            writeVideo(vidfile3,getframe(gcf));
    end
end

close(vidfile3)
  