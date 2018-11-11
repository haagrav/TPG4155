function [gradient, err] = calculateGradient(dims,source,m,trueRec) 
% Inputs -------------------------------------
%
% dims        - Struct containing all relevant dimensions
% m           - Model, velocity field
% source      - Ricker wave source
% trueRec     - True recordings from the real formation
% 
% Outputs -------------------------------------
% 
% gradient    - The gradient of the model
% err         - The error accummilated for all the sources

%% Pre allocations
    u_xy = zeros(dims.nt,length(dims.recPos),'single');    % Recording
    gradient  = zeros(dims.ny,dims.nx,'single');           % Gradient
    u_xy = zeros(dims.my,dims.mx,dims.nt,'single');        % Forwards field 
    u_xy_dagger = zeros(dims.my,dims.mx,dims.nt,'single'); % Backwards field
    err = 0;
    phiComp = 0; % Pre allocation of a comparison Phi. Used in order to print the biggest phi later
    
%% Solving the wave equation    
    for s = 1:dims.ds:length(dims.srcPos)
        
        %% Run forward simulation on background model
        srcPos = dims.srcPos(s);    % Setting the source position to a single source
        u = zeros(dims.ny,dims.nx); % Pre allocating the u, u_1 & u_0 for
        u_0 = u;                    % solving the  wave equation for forward field
        u_1 = u;
        
        for t = 1:dims.nt
            
            % Solve wave equation
            u_1 = waveEqSolver(dims,u,u_0,u_1,m,source,srcPos,t);
            u_0 = u;
            u = u_1;
            
            % Record traces
            u_xt(t,:) = u(dims.recPos);
            
            % Save forward field for use in correlation
            u_xy(:,:,t) = u(dims.modely,dims.modelx);
            
        end % Time loop
        
        %% Calculate difference and error
        
            chi = u_xt - trueRec(:,:,s);
            u_chi = flipud(chi);
            phi = norm(u_chi);
            err = err + sum(abs(u_chi(:)));
            
            %Making sure only the largest phi gets saved
            if phi > phiComp
                phiComp = phi;
            end
        
        %% Run adjoint simulation
        u(:) = 0.0; % Resetting the value of u for the backwards field calculation
        u_0 = u;
        u_1 = u;
        
        chiPos = dims.srcPos; % Setting the source to all source positions
        
        for t = 1:dims.nt
            
            % Backwards propagation. Note that the source is now chi
            u_d = waveEqSolver(dims,u,u_0,u_1,m,u_chi,chiPos,t);
            u_0 = u;
            u = u_d;
            
            % Save forward field for use in correlation
            u_xy_dagger(:,:,dims.nt-t+1) = u(dims.modely,dims.modelx);
            
        end % Time loop
        
        %% Gradient
        for t = 2:dims.nt-1
           
            %Forward field finite difference
            u_xy_M0 = u_xy(:,:,t-1);
            u_xy_M1 = u_xy(:,:,t+1);
            
            u_xy_Dt = (u_xy_M1 - u_xy_M0)./(2*dims.dt);
            
            %Adjoint field finite difference
            u_xy_d_M0 = u_xy_dagger(:,:,t-1);
            u_xy_d_M1 = u_xy_dagger(:,:,t+1);
            
            u_xy_d_Dt = (u_xy_d_M1 - u_xy_d_M0)./(2*dims.dt);
            
            %delta M
            
            gradient(dims.modely,dims.modelx) = gradient(dims.modely,dims.modelx) ...
                                                +( u_xy_d_Dt .* u_xy_Dt);
            
            
            
        end % Time loop
    end %Source loop
    fprintf('The biggest phi is %3.5f\n',phiComp);
end

