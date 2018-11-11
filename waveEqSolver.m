function [u_1] = waveEqSolver(dims,u,u_0,u_1,m,source,srcPos,t)
%  Inputs ------------------------------------
%  
%  dims         - struct containing all relevant dimensions
%  u            - displacement field
%  u_0          - previous displacement field
%  u_1          - the next displacement field
%  m            - model, veolcity field
%  source       - Ricker wave source
%  srcPos       - positions of the source
%  t            - time
%  
%  Outputs ------------------------------------
%  
%  u_1          - the next displacement field
    
    %% Calculation r-values
    r_x = (m.*m).*(dims.dt*dims.dt)/(dims.dx*dims.dx);
    r_y = (m.*m).*(dims.dt*dims.dt)/(dims.dy*dims.dy);
    
    %% Setting the source    
        u(srcPos) = u(srcPos) + source(t,:);
    
    %% Finite difference implementation
        for i=4:dims.ny-3
         for j=4:dims.nx-3
             u_1(i,j) = r_y(i,j).*((1/90)*(u(i-3,j) + u(i+3,j)) ...
                            -(3/20)*(u(i-2,j) + u(i+2,j)) ...
                            + (3/2)*(u(i-1,j) + u(i+1,j)) ...
                            - (49/18)*u(i,j)) ...
                    +r_x(i,j).*((1/90)*(u(i,j-3) + u(i,j+3))...
                            - (3/20)*(u(i,j-2)+ u(i,j+2))...
                            + (3/2)*(u(i,j-1) + u(i,j+1))...
                            - (49/18)*u(i,j)) ...
                            + 2*u(i,j) - u_0(i,j);             
       
         end
       end   
end