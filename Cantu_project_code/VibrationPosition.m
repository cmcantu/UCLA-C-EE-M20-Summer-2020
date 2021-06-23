function[x] = VibrationPosition(x0,m,k,c,f,dt,method)
%function takes in listed variables and approximates position and velocity of
%mass at next time step using designated Runge-Kutta (RK) method
%x0: current state [x0, v0]
%m: mass
%k: stiffness
%c: damping constant
%f: force function 1x3 vector w/ corresponding foreces based on RK1/RK2/RK4
%dt: time step
%method: numerical integration RK1 (Forward Euler), RK2, or RK4

%update variable
x_k = x0(1);
v_k = x0(2);

%calculate spring-damping coefficients
w_n = sqrt(k/m);
xi = c/(2*m*w_n);

%solving with a method
switch method
   %RK1/Forward Euler
    case(1)
       %acceleration (2nd der of x/1st der of v)
        dvdt = -2*xi*w_n*v_k - w_n^2*x_k + f(1);
        
       %solve
        x_(k+1) = x_k + dt*v_k;
        v_(k+1) = v_k + dt*dvdt;
       %store in output vector
        x = [x_(k+1) v_(k+1)];
        
   %RK2
    case(2)
       %acceleration (1st der of v)
        dvdt = -2*xi*w_n*v_k - w_n^2*x_k + f(1);
        
       %constants calculations
        cx1 = dt*v_k;
        cv1 = dt*dvdt;

        cx2 = dt*(v_k+0.5*cv1);
        cv2 = dt*(-2*xi*w_n*(v_k+0.5*cv1) - w_n^2*(x_k+0.5*cx1) + f(2));
       
       %solve
        x_(k+1) = x_k + cx2;
        v_(k+1) = v_k + cv2;
       %store in output vector
        x = [x_(k+1) v_(k+1)];
        
   %RK4   
    case(4)
        %acceleration (1st der of v)
        dvdt = -2*xi*w_n*v_k - w_n^2*x_k + f(1);
        
       %constants calculations
        cx1 = dt*v_k;
        cv1 = dt*dvdt;

        cx2 = dt*(v_k + 0.5*cv1);
        cv2 = dt*(-2*xi*w_n*(v_k + 0.5*cv1) - w_n^2*(x_k + 0.5*cx1) + f(2));
       
        cx3 = dt*(v_k + 0.5*cv2);
        cv3 = dt*(-2*xi*w_n*(v_k+0.5*cv2) - w_n^2*(x_k + 0.5*cx2) + f(2));
        
        cx4 = dt*(v_k + cv3);
        cv4 = dt*(-2*xi*w_n*(v_k + cv3) - w_n^2*(x_k + cx3) + f(3));
        
       %solve
        x_(k+1) = x_k + cx1/6 + cx2/3 + cx3/3 + cx4/6;
        v_(k+1) = v_k + cv1/6 + cv2/3 + cv3/3 + cv4/6;
       %store in output vector
        x = [x_(k+1) v_(k+1)];
        
end
end

