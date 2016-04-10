%_____________________________________________________________________________
%
% Implementation of FDTD in cylindrical coordinate system
% Kapil Saraswat, 13104176
% Electrical Engineering
% Indian Institute of Technology Kanpur
%_____________________________________________________________________________
%=============================================================================

% Program starts from here
%_____________________________________________________________________________

clc;
clear all;
warning off;

% Total no of time steps
%_____________________________________________________________________________
time_tot=200;

% Free space materical charactristics
%_____________________________________________________________________________
epsilon0=1;
mu0=1;
c=1;

% grid step lengths
%_____________________________________________________________________________
deltap=1.8;
deltar=1;
deltat=0.5

% Maximum radius and theta variation
%_____________________________________________________________________________
rmax=30;
tmax=360;

% Variable initialization for source
%_____________________________________________________________________________
T=0;
t0=40.0;
spread=15.0;

% Theta and r stepping 
%_____________________________________________________________________________
rdim=rmax/deltar; 
tdim=tmax/deltap;

% E and H field initialization
%_____________________________________________________________________________
Ez=zeros(rdim,tdim);
Ht=zeros(rdim,tdim);
Hr=zeros(rdim,tdim);

% Initialization of permittivity and permeability matrices
%_____________________________________________________________________________
epsilon=epsilon0*ones(rdim,tdim);
mu=mu0*ones(rdim,tdim);

%_____________________________________________________________________________
% Main FDTD starts from here
%_____________________________________________________________________________

for n=1:deltat:time_tot
    
%_____________________________________________________________________________
% update loops for Ht and Hr fields begin
    for i=1:1:rdim-1
        for j=1:1:tdim-1
            Ht(i,j)=Ht(i,j)+(deltat/(deltar*i*deltap*mu(i,j)))*(Ez(i+1,j)-Ez(i,j));
            Hr(i,j)=Hr(i,j)-(deltat/(deltar*mu(i,j)))*(Ez(i,j+1)-Ez(i,j));
        end
    end
% update loops for Ht and Hr fields end

%_____________________________________________________________________________
% update loop for Ez field begins
    for i=2:1:rdim
        for j=2:1:tdim
            Ez(i,j)=Ez(i,j)+(deltat/(deltar*i*epsilon(i,j)))*((deltar*(i+1/2)/deltar)...
                *Ht(i,j)-(deltar*(i-1/2)/deltar)*Ht(i-1,j)-(1/deltap)*Hr(i,j)+(1/deltap)*Hr(i,j-1));
        end
    end
% update loop for Ez field ends

%_____________________________________________________________________________
% Source conditions at r=0
    
    % pulse=10*sin(pi*0.05*T)           % Sinusoidal Continuous sorce
    pulse=exp(-5.*((t0-T)./spread)^2);  % Gaussian source
    % Update Ez (axial component) for the desired source
    Ez(1,:)=Ez(1,:)+ pulse;
    T=T+1;
  
%_____________________________________________________________________________    
%  Plot in Rectangular form
    %  imagesc(Ez',[-1,1]);colorbar;
    %  title(['Ez plot for 2D polar FDTD (TM)']);
    %  xlabel('r (in spacesteps)');
    %  ylabel('Theta (in Degree)');
    %  getframe;
    
    th = (1:deltap:tmax)*pi/180;
    r = 1:deltar:rmax;
    [TH,R] = meshgrid(th,r);
    [X,Y] = pol2cart(TH,R);
    Z = X + 1i*Y;
    axis([-rmax rmax -rmax rmax -20 20])
    axis tight
    surf(X,Y,abs(Ez))
    view(90,90);
    
%_____________________________________________________________________________
% Contour Plot
    
    % hold on
    % contour(X,Y,Ez,30)
    % colorbar;
    % hold off;
    
    getframe;

end   
%_____________________________________________________________________________
% End of program
%_____________________________________________________________________________

% Effect of 1/r
%_____________________________________________________________________________

% for kk=1:rmax
% plot(Ez(:,kk));
% end
%_____________________________________________________________________________