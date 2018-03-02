function tmpNext = heat_diffusion_CN(tmpCurr, bcTyp1, bcVal1, bcTypN, bcValN, c)

nLyr = numel(tmpCurr) - 1;

%Initialize matrix system to be solved:
matCf = zeros([nLyr+1, nLyr+1], 'single');

%Assign a, b, and c:
matCf(1:nLyr+2:end) = 2*c+1; %main diagonal
matCf(nLyr+2:nLyr+2:end) = -c; %upper diagonal
matCf(2:nLyr+2:end) = -c; %lower diagonal

vecKwn = nan([nLyr+1, 1]);
vecKwn(2:end-1) = c*tmpCurr(1:end-2) + (1-2*c)*tmpCurr(2:end-1) + c*tmpCurr(3:end);

%Set BC
if strcmpi(bcTyp1, 'dirichlet') %Constant temperature
    matCf(1,1) = 1;
    matCf(1,2) = 0;
    vecKwn(1) = bcVal1;
elseif strcmpi(bcTyp1, 'neumann') %Constant heat flux
    %matCf(1,1) = 2*c+1; %already assigned
    matCf(1,2) = 0;
    vecKwn(1) = (1-2*c)*tmpCurr(1) + bcVal1*c;
    %BC Val = 2*deltaZ*(qPrev+qCurr)/k
else
    error('heatDiffusionCn:bc1Unknown', ['The boundary condition ' bcTyp1 ' has not been programmed for.']);
end

if strcmpi(bcTypN, 'dirichlet') %Constant temperature
    matCf(end,end) = 1;
    matCf(end,end-1) = 0;
    vecKwn(end) = bcValN;     
elseif strcmpi(bcTypN, 'neumann') %Constant heat flux
    %matCf(end,end) = 2*c+1; %already assigned
    matCf(end,end-1) = 0;
    vecKwn(1) = (1-2*c)*tmpCurr(end) + bcValN*c;
else
    error('heatDiffusionCn:bc1Unknown', ['The boundary condition ' bcTypN ' has not been programmed for.']);
end

%keyboard

%Solve:
tmpNext = linsolve(matCf, vecKwn);



% %FROM http://web.cecs.pdx.edu/~gerry/class/ME448/notes/pdf/CN_slides.pdf
% % --- Coefficients of the tridiagonal system
% b = (-alfa/2/dx^2)*ones(nx,1); % Super diagonal: coefficients of u(i+1)
% c = b; % Subdiagonal: coefficients of u(i-1)
% a = (1/dt)*ones(nx,1) - (b+c); % Main Diagonal: coefficients of u(i)
% at = (1/dt + b + c); % Coefficient of u_i^k on RHS
% a(1) = 1; b(1) = 0; % Fix coefficients of boundary nodes
% a(end) = 1; c(end) = 0;
% [e,f] = tridiagLU(a,b,c); % Save LU factorization
% % --- Assign IC and save BC values in ub. IC creates u vector
% x = linspace(0,L,nx)'; u = sin(pi*x/L); ub = [0 0];
% % --- Loop over time steps
% for k=2:nt
%     % --- Update RHS for all equations, including those on boundary
%     d = - [0; c(2:end-1).*u(1:end-2); 0] ...
%     + [ub(1); at(2:end-1).*u(2:end-1); ub(2)] ...
%     - [0; b(2:end-1).*u(3:end); 0];
%     u = linsolve(A,B); % Solve the system
% end
% 
% 
% %FROM http://www.goddardconsulting.ca/matlab-finite-diff-crank-nicolson.html
% % Form the tridiagonal matrix
% C = -diag(a(3:nLayers),-1) + diag(1-b(2:nLayers)) - diag(c(2:nLayers-1),1);
% [L,U] = lu(C);
% D = diag(a(3:nLayers),-1) + diag(1+b(2:nLayers)) + diag(c(2:nLayers-1),1);
% 
% % Solve at each node
% offset = zeros(size(D,2),1);
% for idx = N:-1:1
%     if length(offset)==1
%         offset = a(2)*(price(1,idx)+price(1,idx+1)) + ...
%             c(end)*(price(end,idx)+price(end,idx+1));
%     else
%         offset(1) = a(2)*(price(1,idx)+price(1,idx+1));
%         offset(end) = c(end)*(price(end,idx)+price(end,idx+1));
%     end
%     price(2:nLayers,idx) = U\(L\(D*price(2:nLayers,idx+1) + offset));
% end
