%% Clearing the work space
close all;
clear all;
clc

%% Setting up dimensions
dims.dy =     10; % [m]
dims.dx = dims.dy;% [m]
dims.dt = 1.0e-3; % [s]

dims.ny = 201; % Cells in y-direction
dims.nx = 301; % Cells in x-direction
dims.nt = 801; % Amount of time steps

%% Model dimensions
dims.modely = 100:150;
dims.modelx = 100:200;
dims.my = length(dims.modely);
dims.mx = length(dims.modelx);

%% Source locations
sx = min(dims.modelx):max(dims.modelx);
sy = min(dims.modely)*ones(1,length(sx));
dims.srcPos = sy + dims.ny*sx;

%% Receiver locations
rx = min(dims.modelx):max(dims.modelx);
ry = min(dims.modely)*ones(1,length(rx));
dims.recPos = ry+dims.ny*rx;

%% Creating background model
bg = zeros(dims.ny,dims.nx,'single');
bg(:) = 2.0e3;         % [m/s] - Background
bg(115:end,:) = 2.3e3; % [m/s] - Layer

%% Begin iteration
m = bg;                 % Starting model
dims.ds = 10;           % Grid point distance between sources
maxIter = 10;           % Maximum number of iterations per frequency
freqs = [4,6,8,10,12];  % Frequencies to use in inversion
errVec = zeros(1,maxIter*length(freqs));
      
%% Full Waveform Inversion
 for f=freqs
 it = 1;
 
 source = rickerWave(f,dims);
 load (['trueRec_',num2str(f),'Hz.mat']); 
 tic;
 
 for i = 1:maxIter

     fprintf('----Iteration %d------\n',it);
     fprintf('The frequency is %d\n',f);
     
     %% Calculating the gradient and error
     fprintf('(Calculating gradient)\n');
     [gradient,err] = calculateGradient(dims,source,m,trueRec);
     fprintf('The error for iteration %d is %3.5f\n',it,err);
     
     %% Tapering the gradient
     gradient = taperGradient(gradient);
     
     %% Calculating the step length
     fprintf('(Calculating step length)\n');
     [stepLength,newErr] = calculateStepLength(dims,gradient,err,source,m,trueRec);
     
     %% Updates
     
     % Updating the model
     m = m + gradient*stepLength;
     errVec(it) = err;
     
     %Plotting the new model
     clf;
     imagesc(m(dims.modely,dims.modelx)); colorbar;
     drawnow;
     shg;
     
     toc;
     it=it+1;
 end
 end       