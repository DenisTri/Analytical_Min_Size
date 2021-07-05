function[R, muEro, muDil, tDil, tEro] = SizeSolution(eta_ero, rSol, rVoid, beta, epsVar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to find the parameters that induces minimum sizes               %
% R is the radius filter                                                   %
% muEro is the erode threshold                                             %
% muDil is the dilated threshold                                           %
% rSol is the radius of the minimum size in the solid phase                %
% rVoid  the radius of the minimum size in th void phase                   %
% beta is the final parameter of the heaviside projection                  %
% epsVar define the cut-off threshold                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta_i = 0.5;                               % Intermediate threshold

%% Solid
%eta_ero = 0.55:0.01:0.95;                    % Eroded threshold
eps = epsVar;

eta = eta_i + 1/beta*atanh(2*eps-1);       % Correct intermediate threshold
eta_e = eta_ero + 1/beta*atanh(2*eps-1);   % Correct eroded threshold

% Compute the different zones possible (relation between eta_i and eta_e)
[Zone1, Zone2, Zone3, Zone4] = SolZone(eta, eta_e);
% In each zone, compute the normalized minimum size
b_Sol = minSol(eta, eta_e, Zone1, Zone2, Zone3, Zone4);

Rfil = 2*rSol./b_Sol; % Compute the filter radius

R = Rfil';
muEro = eta_ero';

%% Void
eps = 1-epsVar;
eta = eta_i + 1/beta*atanh(2*eps-1);
eta_d = zeros(1,length(eta_e));

% Compute the dilated threshold for all zones and radius filter
d1 = (sqrt(2*eta)-rVoid./Rfil).^2;     
d2 = eta - (rVoid./Rfil).^2;           
d3 = (2-rVoid./Rfil-sqrt(2-2*eta)).^2; 
d4 = (1+(1-eta)./(rVoid./Rfil-2)).^2;  

% Find and assing the different zones (relation btw eta_e and eta_i)
[Zone1,~, ~, ~] = VoidZone(eta, d1); eta_d(Zone1) = d1(Zone1);
[~,Zone2, ~, ~] = VoidZone(eta, d2); eta_d(Zone2) = d2(Zone2);
[~, ~,Zone3, ~] = VoidZone(eta, d3); eta_d(Zone3) = d3(Zone3);
[~, ~, ~,Zone4] = VoidZone(eta, d4); eta_d(Zone4) = d4(Zone4);

eta_d = eta_d - 1/beta*atanh(2*eps-1); % Correct the dilated threshold


%% Validity interval
Validity = find(eta_d<0.01);
eta_ero(Validity:end) = [];
eta_d(Validity:end) = [];
eta_e(Validity:end) = [];

b_Sol(Validity:end) = [];
Rfil(Validity:end) = [];
R(Validity:end) = [];
muEro(Validity:end) = [];

muDil = eta_d';
%% t_dil

eps = epsVar;
eta_d2 = eta_d + 1/beta*atanh(2*eps-1);

% Compute the different zones possible (considering eta_d as "intermediate")
[Zone1, Zone2, Zone3, Zone4] = SolZone(eta_d2, eta_e);

% In each zone, compute the minimum size
b_d = minSol(eta_d2, eta_e, Zone1, Zone2, Zone3, Zone4);

tDil = (b_d-b_Sol)/2 .* Rfil;

tDil = tDil';

%% t_ero

eps = 1-epsVar;
eta_e = eta_ero + 1/beta*atanh(2*eps-1);
eta_d3 = eta_d + 1/beta*atanh(2*eps-1);

[Zone1, Zone2, Zone3, Zone4] = VoidZone(eta, eta_d3);
b_V = minVoid(eta, eta_d3, Zone1, Zone2, Zone3, Zone4);

[Zone1, Zone2, Zone3, Zone4] = VoidZone(eta_e, eta_d3);
b_e = minVoid(eta_e, eta_d3, Zone1, Zone2, Zone3, Zone4);

tEro = (b_e-b_V)/2 .* Rfil;

tEro = tEro';
end

function[Zone1, Zone2, Zone3, Zone4] = VoidZone(eta, eta_d)
% Function compute the different zones
% eta is the intermediate threshold
% eta_d is a scalar/vector of dilated threshold

% Zone 1 (b<h=h, R<=h)
condition1 = eta<=1/2;
condition2 = eta>=2*eta_d;
Zone1 = condition1 & condition2;

% Zone 2 (b<h=h, R>h)
condition3 = eta<2*eta_d;
condition4 = eta<=2*eta_d + 1 - 2*sqrt(eta_d);
Zone2 = condition3 & condition4;

% Zone 3 (b>h=h, R<h)
condition5 = eta>1/2;
condition6 = eta>=4*sqrt(eta_d) - 2*eta_d - 1;
Zone3 = condition5 & condition6;

% Zone 4 (b>h, R>=h)
condition7 = eta>2*eta_d + 1 - 2*sqrt(eta_d);
condition8 = eta<4*sqrt(eta_d) - 2*eta_d - 1;
Zone4 = condition7 & condition8;
end

function[Zone1, Zone2, Zone3, Zone4] = SolZone(eta, eta_e)
% Function compute the different zones
% eta is the intermediate threshold
% eta_e is a scalar/vector of eroded threshold

% Zone 1 (b<h=h, R<=h)
condition1 = eta>=1/2;
condition2 = eta<=2*eta_e-1;
Zone1 = condition1 & condition2;

% Zone 2 (b<h=h, R>h)
condition3 = eta>2*eta_e-1;
condition4 = eta>=2*eta_e-2+2*sqrt(1-eta_e);
Zone2 = condition3 & condition4;

% Zone 3 (b>h=h, R<h)
condition5 = eta<1/2;
condition6 = eta<4-4*sqrt(1-eta_e) - 2*eta_e;
Zone3 = condition5 & condition6;

% Zone 4 (b>h, R>=h)
condition7 = eta>=4-4*sqrt(1-eta_e) - 2*eta_e;
condition8 = eta<2*eta_e-2+2*sqrt(1-eta_e);
Zone4 = condition7 & condition8;
end

function[rMinVoid] = minVoid(eta, eta_d, Zone1, Zone2, Zone3, Zone4)
% Function compute minimum size in the void phase
% eta is the corrected intermediate threshold
% eta_d is a scalar/vector of corrected dilated threshold

if length(eta)==1
    eta = ones(size(eta_d))*eta;
end

b_R = zeros(1,length(eta_d));

% In each zone, compute the minimum size
b_R(Zone1) = 2*sqrt(2*eta(Zone1)) - 2*sqrt(eta_d(Zone1));
b_R(Zone2) = 2*sqrt(eta(Zone2)-eta_d(Zone2));
b_R(Zone3) = 4 - 2*sqrt(eta_d(Zone3)) - 2*sqrt(2-2*eta(Zone3));
b_R(Zone4) = 2 - (1-eta(Zone4))./(1-sqrt(eta_d(Zone4)));

rMinVoid = b_R;
end

function[rMinSol] = minSol(eta, eta_e, Zone1, Zone2, Zone3, Zone4)
% Function compute minimum size in the void phase
% eta is the corrected intermediate threshold
% eta_e is a scalar/vector of corrected eroded threshold

if length(eta)==1
    eta = ones(size(eta_e))*eta;
end

b_R = zeros(1,length(eta_e));

% In each zone, compute the minimum size
b_R(Zone1) = 2*sqrt(2-2*eta(Zone1)) - 2*sqrt(1-eta_e(Zone1));
b_R(Zone2) = 2*sqrt(eta_e(Zone2)-eta(Zone2));
b_R(Zone3) = 4 - 2*sqrt(1-eta_e(Zone3)) - 2*sqrt(2*eta(Zone3));
b_R(Zone4) = 2 - eta(Zone4)./(1-sqrt(1-eta_e(Zone4)));

rMinSol = b_R;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Fernandez, D. Trillet                 %
% Dept. of Aerospace and Mechanics, University of Liège                    %
% 4000 Liège (BEL), July 2020                                              % 
% Please send your comments to: dtrillet@uliege.be or efsanchez@uliege.be  %
%                                                                          %
% Theoretical details are discussed in the paper:                          %
% Note on the minimum length scale and its defining parameters, 2020       %
% Trillet D., Fernandez E., Duysinx, P.                                    %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
