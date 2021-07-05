function [] = MinLengthScales_Analytical_and_Numerical(rfil)
  
% This code computes numerically the minimum size in the solid and void phases
% The offset distances are computed as well. All the datas are saved in
% MinVoid.m and  MinSolid.m. 
% The only argument is the filter radius. See the function "parameters.m" to 
% change the beta parameter or other.
% The graphs are related to the numerical part of the graphs 6.a) to 6.d) 
% presented in the paper :
% Note on the minimum length scale and its defining parameters, 2020
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

 if ~exist('rfil','var')
     % If the filter radius is not declare, set it to 40
     rfil = 40;
 end

%% ========================================================================
%========================== NUMERICAL ===================================== 
%% ========================================================================
%================ MININIMUM SIZE IN THE SOLID PHASE =======================

muint = 0.5;             % threshold for the intermediate field
EpsiSol = 0.5;
if exist('MinSolid.MAT', 'file')==0 % computing the data for MinSolid
  fprintf(1,'\n --------- Computing MinSolid ---------');
  [rfil,beta,N,MU,rInt,tdil,H,nMU] = Parameters(rfil);    
  for c=1:nMU
    mu   = MU(c)+muint;            % threshold of the eroded field
    % This corr value can be used for speeding up the While loop.
    corr = 1;
    epsi = EpsiSol;                % Considered as Solid element
    x    = zeros(N,1);             % Field to be filtered
    x(N/2-round(corr/2):N/2-round(corr/2)+corr) = 1;      % Filling x with solid elements
    xfil = H*x;                    % Filtered field
    Ero  = Heaviside(xfil,mu,beta);% Eroded field
    [iEro,~] = find(Ero>epsi);     % Check the amount of solid elements
    if isempty(iEro)==1            
      while isempty(iEro)==1       % Enter if No solid elements
        corr = corr+1;             % Increase the amount of solid elements
        x    = zeros(N,1);         % The field
        x(N/2-round(corr/2):N/2-round(corr/2)+corr) = 1;        % The increased amount of solid
        xfil = H*x;                      % The filtered field
        Ero   = Heaviside(xfil,mu,beta); % The eroded field
        % To check the amount of solid.
        [iEro,~] = find(Ero>=epsi);
      end
    end
    Int      = Heaviside(xfil,muint,beta); % The intermediate field
    % The solid inside Int.
    [iInt,~] = find(Int>=epsi);   
    rInt(c,1)= (iInt(end)-iInt(1)+1)/rfil/2; % The size of Int. (RminSolid)
    % Loop over muDil ------------------------------------------------||||
    for cc = 1:nMU
      mudil = MU(cc);                         % The dilation threshold
      Dil   = Heaviside(xfil,mudil,beta);     % The dilated fiel
      % The Solid portion
      [iDil,~] = find(Dil>=epsi); 
      rDil  = (iDil(end)-iDil(1)+1)/rfil/2;     % The size of Dil (RminSolidDil)
      tdil(cc,c) = abs(rDil-rInt(c,1))/rInt(c,1);% The offset distance (tdil)
    end % ------------------------------------------------------------||||
    if rem(c,2)==0,fprintf(1,'\n Progress(Solid): %3.0f%% \n',c/nMU*100);end  
  end  
  MinSolid = {MU,rInt,tdil}; save('MinSolid','MinSolid'); % Save Info.
end

%% ==========================================================================
%================= MININIMUM SIZE IN THE VOID PHASE   =======================
if exist('MinVoid.MAT', 'file')==0         % computing the data for MinVoid
  fprintf(1,'\n --------- Computing MinVoid ---------\n');
  [rfil,beta,N,MU,rInt,tero,H,nMU] = Parameters(rfil);
  for c=1:nMU
    mu   = MU(c);              % threshold of the dilated field
    % This corr value can be used for speeding up the While loop.
    corr = 1;
    epsi = 1-EpsiSol;               % Considered as Void element
    x    = ones(N,1);          % The field to be filtered 
    x(N/2-round(corr/2):N/2-round(corr/2)+corr) = 0;  % Filling the field with void elements 
    xfil = H*x;                % Filtered field 
    Dil  = Heaviside(xfil,mu,beta); % Dilated field 
    [iDil,~] = find(Dil<=epsi); % Dilated field 
    if isempty(iDil)==1        % No void elements in Dil
      while isempty(iDil)==1   % then, increase the void in x
        corr = corr+1;         % 1 more void element (radius)
        x    = ones(N,1);      % the field
        x(N/2-round(corr/2):N/2-round(corr/2)+corr) = 0; % the field with void elements
        xfil = H*x;            % the filtered field
        Dil   = Heaviside(xfil,mu,beta); % The dilated field
       % Finding void element
        [iDil,~] = find(Dil<=epsi); % Dilated field 
      end
    end
      Int      = Heaviside(xfil,muint,beta); % The intermediate field
      % Finding void elements
      [iInt,~] = find(Int<=epsi);  
      rInt(c,1)= (iInt(end)-iInt(1)+1)/rfil/2; % The minimum size (void)
      % Loop over muEro -----------------------------------------------||||
      for cc = 1:nMU                         
        mu       = MU(cc)+muint;                     % The ersion threshold
        Ero      = Heaviside(xfil,mu,beta);         % The eroded field
        %[iEro,~] = find(xfil<mu);                  % The void elements
        [iEro,~] = find(Ero<=epsi);  
        rEro     = (iEro(end)-iEro(1)+1)/rfil/2;      % the minimum size (void)
        tero(cc,c) = abs(rEro-rInt(c,1))./rInt(c,1);% The offset distance tero
      end % -----------------------------------------------------------||||
      if rem(c,2)==0,fprintf(1,' Progress(Void): %3.0f%%\n',c/nMU*100);end  
  end
  MinVoid = {MU,rInt,tero}; save('MinVoid','MinVoid'); % Save info
end

%% Plot the figures relative to the data just computed
load('MinSolid');rIntSol=MinSolid{2}; muEro=MinSolid{1}+0.5; tdil=MinSolid{3};
load('MinVoid'); rIntVoi=MinVoid{2};  muDil=MinVoid{1};      tero=MinVoid{3};  

% rMinSol/rMinVoid vs muEro
figure;
hold on;
for c= 1:size(rIntVoi,1)
    F1a=plot(muEro,rIntSol./rIntVoi(c,1),'b-x','linewidth',0.5,'MarkerSize',5,'MarkerEdgeColor','r'); hold on;
    text(0.96,min(rIntSol(end)/rIntVoi(c)-0.01,3.45),num2str(muDil(c)),'color','b','fontsize',12);
end; axis([0.5 1.0 0 3.5]); set(gca,'fontsize',16);
text(0.95,0.40,'$\eta_\mathrm{dil}$',...
'interpreter','latex','color','b','fontsize',17); annotation('arrow',... 
[0.86 0.86],[0.28 0.13+rIntSol(end)/rIntVoi(1)/5],'color','b');
Labels('$\eta_\mathrm{ero}$','$\frac{r_\mathrm{min.Solid}^\mathrm{int}}{r_\mathrm{min,Void}^\mathrm{int}}$',18,22);

% rMinSol vs muEro
figure; hold on;
F2=plot(muEro, 2*rIntSol,'-x','linewidth',1.2,'MarkerSize',5,'MarkerEdgeColor','r');

set(gca,'xtick',[0.5:0.1:1.])
set(gca,'ytick',[0.2:0.2:1.8])

axis([0.5 1.0 0.2 1.8]); set(gca,'fontsize',18);
Labels('$\eta_\mathrm{ero}$','$\frac{2r_\mathrm{min,Sol}^\mathrm{int}}{r_\mathrm{fil}}$', 18, 22);

% tDil vs muEro (not accurate for rfil < 30)
figure
hold on

for c = 1:size(tdil,1)
   F3a=plot(muEro',tdil(c,:),'b-x','linewidth',0.5,'MarkerSize',5,'MarkerEdgeColor','r'); hold on;
   text(0.96,tdil(c,end)+0.01,num2str(muDil(c)),'color','b','fontsize',12);
end; axis([0.5 1.0 0 2]); set(gca,'fontsize',16);
text(0.955,1.20,'$\eta_\mathrm{dil}$',...
'interpreter','latex','color','b','fontsize',18); annotation('arrow',... 
[0.86 0.86],[0.59 0.54],'color','b');
Labels('$\eta_\mathrm{ero}$','$\frac{t_\mathrm{dil}}{r_\mathrm{min,Solid}^\mathrm{int}}$',18,24)

% tEro vs muEro (not accurate for rfil < 30)
figure
hold on

for c = 1:size(tdil,1) 
 F4a=plot(muEro',tero(:,c),'b-x','linewidth',0.5,'MarkerSize',5,'MarkerEdgeColor','r'); hold on;
   if tero(end,c) <3
      text(0.96,tero(end,c)+0.01,num2str(muDil(c)),'color','b','fontsize',12);
   end
end; axis([0.5 1 0 3]); set(gca,'fontsize',16);
text(0.95,0.35,'$\eta_\mathrm{dil}$',...
'interpreter','latex','color','b','fontsize',18); annotation('arrow',... 
[0.86 0.86],[0.28 0.13+tero(end,1)/4],'color','b');

Labels('$\eta_\mathrm{ero}$','$\frac{t_\mathrm{ero}}{r_\mathrm{min,Void}^\mathrm{int}}$',18,24)


  function H = Heaviside(x,mu,beta)
  % Function to compude the Heaviside projection of a field
  H=(tanh(beta*mu)+tanh(beta*(x-mu)))/(tanh(beta*mu)+tanh(beta*(1-mu)));

  function [rfil,beta,N,MU,rInt,t,H,nMU] = Parameters(rfil)
  % Function to choose the parameters of the reference optimization   
  beta    = 258;                 % For heaviside function
  N       = 4*rfil;              % Size of Design domain
  coord   = (1:N)';              % Coordinates to create de Filter Matrix
  MU      = [(0.05:0.05:0.45)']; % Reduce step size to improve Resolution
  nMU     = length(MU);          % Amount of thresholds
  rInt    = zeros(nMU,1);        % To store minSize of Intermediate design
  t       = zeros(nMU,nMU);      % To store the offset distances 
  IND 	  = cell(N,1);           % ------------------------- FILTER ----%
  for el = 1:N                                                          %
    dist  = abs(coord(el)-coord);                                       %
    [i,j] = find(dist<=rfil);                                           %
    IND{el,1} = [el+j-1, i, (1-dist(i)/rfil)];                          %
  end                                                                   %
  IND = cell2mat(IND); H   = sparse(IND(:,1),IND(:,2),IND(:,3));        %
  sumH= H*ones(N,1); H = spdiags(1./sumH,0,N,N)*H; %------ end Filter --% 

  function [] = Labels(X,Y, x, y)
  xlabel(X,'interpreter','latex','fontsize',x); %grid minor; %23
  ylabel(Y,'interpreter','latex','fontsize',y); %32

