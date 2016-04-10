% Assignment 4; EE640
% File Name: Assignment4.m
% Kapil Saraswat,13104176
% Department of Electrical Engineering
% Indian Institute of Technology Kanpur

clc;
close all;
clear all;
warning off;
tic;

dx=150;                 % number of nodes in x direction
dy=150;                 % number of nodes in y direction
xc=fix(dx/2);           % Source location at middle X-direction
yc=fix(dy/2);           % Source location at middle Y-direction

nsteps=250;             % No. of time steps=250
t0=40.0;
spread=15.0;
npml=20;                % no. of PML nodes  in left, right, top and upper

% Initialization of E,H ans D matrix
ez=zeros(dx,dy);
dz=zeros(dx,dy);
hy=zeros(dx,dy);
ihx=zeros(dx,dy);
ihy=zeros(dx,dy);
hx=zeros(dx,dy);
ga=ones(dx,dy);

% See comment 1 for the initialization of PML
%Initialization of PML parameter in X direction
gi2=ones(dx);
gi3=ones(dx);
fi1=zeros(dx);
fi2=ones(dx);
fi3=ones(dx);

%Initialization of PML parameter in Y direction
gj2=ones(dy);
gj3=ones(dy);
fj1=zeros(dy);
fj2=ones(dy);
fj3=ones(dy);
    
% Calculation of PML parameters for left and right boundaries
for i=1:npml
    xnum=npml-i;
    % See comment2 for the calculation of xn
    xn=0.33*((xnum/npml)^3);        % see equation 3.20
    gj2(i)=1.0/(1.0+xn);            % see equation 3.21b
    gi2(dx-1-i)=1.0/(1.0+xn);       % see equation 3.21b
    gi3(i)=(1.0-xn)/(1.0+xn);       % see equation 3.21c
    gi3(dx-i-1)=(1.0-xn)/(1.0+xn);  % see equation 3.21c
    fi1(i)=xn;                      % see equation 3.21a
    fi1(dx-2-i)=xn;                 % see equation 3.21a
    fi2(i)=1.0/(1.0+xn);            % see equation 3.21b
    fi2(dx-2-i)=1.0/(1.0+xn);       % see equation 3.21b
    fi3(i)=(1.0-xn)/(1.0+xn);       % see equation 3.21c
    fi3(dx-2-i)=(1.0-xn)/(1.0+xn);  % see equation 3.21c
end
% Calculation of PML parameters for top and bottom boundaries
for j=1:npml
    xnum=npml-j;
    xn=0.33*((xnum/npml)^3);        % see equation 3.20
    gj2(j)=1.0/(1.0+xn);            % see equation 3.21b
    gi2(dy-1-j)=1.0/(1.0+xn);       % see equation 3.21b
    gj3(j)=(1.0-xn)/(1.0+xn);       % see equation 3.21c
    gj3(dy-j-1)=(1.0-xn)/(1.0+xn);  % see equation 3.21c
    fj1(j)=xn;                      % see equation 3.21a
    fj1(dy-2-j)=xn;                 % see equation 3.21a
    fj2(j)=1.0/(1.0+xn);            % see equation 3.21b
    fj2(dy-2-j)=1.0/(1.0+xn);       % see equation 3.21b
    fj3(j)=(1.0-xn)/(1.0+xn);       % see equation 3.21c
    fj3(dy-2-j)=(1.0-xn)/(1.0+xn);  % see equation 3.21c
end
T=0;
% Implementation of FDTD Algorithm
for n=1:nsteps
    T=T+1;
    for j=2:dy
        for i=2:dx
            % For 2-D; see equation 3.14 
            dz(i,j)=gi3(i)*gj3(j)*dz(i,j)+gi2(i)*gj2(j)*0.5*(hy(i,j)-hy(i-1,j)...
                -hx(i,j)+hx(i,j-1));
        end
    end
    % Source
    %pulse=10*sin(pi*0.05*T)         % Sinusoidal Continuous sorce
    pulse=exp(-5.*((t0-T)./spread)^2);  % Gaussian sorce
    dz(xc,yc)=dz(xc,yc)+pulse;
    for j=2:dy
        for i=2:dx
            ez(i,j)=ga(i,j).*dz(i,j);   % see equation 3.8b
        end
    end
    for j=1:dy-1
        for i=1:dx-1
            % see equation 3.18a; X-direction
            curl_e=ez(i,j)-ez(i,j+1);
            % see equation 3.18b; X-direction
            ihx(i,j)=ihx(i,j)+fi1(i)*curl_e;
            % see equation 3.18c; X-direction
            hx(i,j)=fj3(j)*hx(i,j)+fj2(j)*0.5*(curl_e+ihx(i,j));
        end
    end
       
     for j=1:dy-1
        for i=1:dx-1
            % see equation 3.18a; Y-direction
            curl_e=ez(i+1,j)-ez(i,j);
            % see equation 3.18b; Y-direction
            ihy(i,j)=ihy(i,j)+fj1(j)*curl_e;
            % see equation 3.18c; Y-direction
            hy(i,j)=fi3(i)*hy(i,j)+fi2(i)*0.5*(curl_e+ihy(i,j));
        end
     end
 colormap(jet);
 %colormap(gray);
 pppp=mesh(ez);
 view(90,90);
 title(['Time = ',num2str(n)]);
 xlabel(' length')
 ylabel(' length')
 axis([0 150 0 150 -2 4])
 axis square
 pause(0.002);
 %eval(['print -djpeg print_npml' num2str(n) '.jpeg']);
end