function [x,y] = RandWalk_2D(x0, y0, BC)
%%% function simulates random walker moving/staying put in 2D grid
%takes in current x & y positions and boundary condition vector(BC)
%walker can move up, down, left, right, or stay put & cannot move past BC
%BC: [up,down,left,right]

%generate random move option
r = rand;

%move up
if r<=0.2
    x=x0; y=y0+1;
    %check upper boundary
    if y>=BC(1)
        y=BC(1);
    end
    
%move down
elseif 0.2<r && r<=0.4
    x=x0; y=y0-1;
    %check lower boundary
    if y<=BC(2)
        y=BC(2);
    end
    
%move left
elseif 0.4<r && r<=0.6
    x=x0-1; y=y0;
    %check left boundary
    if x<=BC(3)
        x=BC(3);
    end
    
%move right
elseif 0.6<r && r<=0.8
    x=x0+1; y=y0;
    %check right boundary
    if x>=(4)
        x=BC(4);
    end
%stay put
elseif 0.8<r
    x=x0; y=y0;
end
end

    