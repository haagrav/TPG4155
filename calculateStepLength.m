function [stepLength,newErr] = calculateStepLength(dims,gradient,oldErr,source,m,trueRec)
%  Inputs-------------------------------------    
%  
% dims            - struct containing all relevant dimensions
% oldErr          - the error acquired from the calculateGradient function
% source          - Ricker wave source
% m               - model, velocity field
% trueRec         - True recordings from the real formation
%
% Outputs ------------------------------------
% 
% stepLenght     - step length
% err            - the error accumilated

%% Pre allocating
    u_reca = zeros(dims.nt,length(dims.recPos),'single'); % Recordings for alpha
    u_recb = zeros(dims.nt,length(dims.recPos),'single'); % Recordings for beta
    stepLength =500;
    newErr = inf;
    
    a = 0;              % Starting point, arbitrairily chosen
    b= 500;             % End point, arbitrairily chosen
    r = (-1+sqrt(5))/2; % The golden ratio
    GSi = 1;
    
%% Golden Search Method        
while (newErr > oldErr)
    
    fprintf('-----Golden search iteration %d -----\n',GSi);
    
    %% Pre allocations
        error_a= 0;
        error_b=0;
        newErr = 0;
        
        phi_a = 0;
        phi_A = 0;
        
        phi_b = 0;
        phi_B = 0;
        
        %% Getting the alpha and beta values
        alpha = a + (1-r)*(b-a);
        beta = b - (1-r)*(b-a);
        fprintf('Alpha = %3.5f & Beta = %3.5f\n',alpha,beta);
               
        %% Setting the test velocityfields
        ma = m + gradient*alpha;
        mb = m + gradient*beta;
        
        %% Solving the wave equation for ma and mb 
    for s = 1:dims.ds:length(dims.srcPos)    
        
        srcPos = dims.srcPos(s);
        
        % Pre allocating the displacement fields
        ua = zeros(dims.ny,dims.nx);
        u_1a = ua;
        u_0a = ua;
        ub = zeros(dims.ny,dims.nx);
        u_1b = ub;
        u_0b = ub;
        
        for t = 1:dims.nt
            
         u_1a = waveEqSolver(dims,ua,u_0a,u_1a,ma,source,srcPos,t);
         u_0a = ua;
         ua = u_1a;
         
         % Recording traces
         u_reca(t,:) = ua(dims.recPos);
         
         
         u_1b = waveEqSolver(dims,ub,u_0b,u_1b,mb,source,srcPos,t);
         u_0b = ub;
         ub = u_1b;
         
         % Recording traces
         u_recb(t,:) = ub(dims.recPos);
         
        end % Time loop
        
        %% Getting the error and chi values
        chi_a = u_reca - trueRec(:,:,s);
        chi_b = u_recb - trueRec(:,:,s);
        
        u_chi_a = flipud(chi_a);
        u_chi_b = flipud(chi_b);
        
        error_a = error_a + sum(abs(u_chi_a(:)));
        error_b = error_b + sum(abs(u_chi_b(:)));
        
        phi_a = norm(u_chi_a);
        phi_b = norm(u_chi_b);
        
        if phi_a> phi_A
            phi_A = phi_a;
        end
        if phi_b> phi_B
            phi_B = phi_b;
        end  
        
    end % End Source
    
    %% Updating the Golden Search Method
if phi_A < phi_B
    b = beta;
    newErr = error_a;
    fprintf('Alpha = %3.5f is the best steplength\n',alpha);
    fprintf('-------------------------------------\n');
    stepLength = alpha;
else
    a = alpha;
    newErr = error_b;
    fprintf('Beta = %3.5f is the best steplength\n',beta);
    fprintf('-------------------------------------\n');
    stepLength = beta;
end
    GSi = GSi + 1;
    
    % End criterium to the while loop
    if GSi ==15
        newErr = 0;
    end
end % End while
end
