

clear
close all

% base = [-1 -1 1 1; 0 -.2 -.2 0]/10;
base = [-1 -1 1 1; .2 0 0 0.2]/4;

%{
C = interpft(randn(10,3), 501)'; %Generate a smooth periodic path
C = 1*C./max(abs(C(:))); %normalize it
%}

% Create a test curve
s = linspace(0, pi, 101);
a = 0.5;
gam = (s+a).*exp(1i*s);

ft = @(z) z.^2;
f = ft(gam);

C = [gam(:).'; f(:).'];
% C = 1*C./max(abs(C(:))); %normalize it

% Use 'real', 'imag', 'abs', etc.
reimname = 'real';
reim = str2func(reimname);

xx = linspace(-4, 4, 50);
yy = xx;
[XX,YY] = meshgrid(xx,yy);

ZZ = reim(ft(XX + 1i*YY));
figure
% surf(XX, YY, ZZ-0.01, 'FaceColor', 'b', 'FaceAlpha', 0.5); hold on;
% hold on

%% Calculate derivatives and allocate matrices
npt = size(base,2);
base = [base; zeros(1,npt)];

if size(C,2) >= 3 %Use a 2nd order approximation for the derivatives of the trajectory
    dC = [C(:,1:3)*[-3; 4; -1]/2  [C(:,3:end) - C(:,1:end-2)]/2 C(:,end-2:end)*[1; -4; 3]/2];   
else
    dC = C(:,[2 2]) - C(:,[1 1]);
end

dC0 = find(sum(abs(dC),1) == 0,1);    %Check for stagnation points
if ~isempty(dC0)
    warning('Removing stagnation points found in trajectory');
    dCgood = find(sum(abs(dC),1) ~= 0);
    C = C(:,dCgood); %
    dC = dC(:,dCgood); %
end

K = size(dC,2);
SUR = nan(3,npt,K);

%% Generate and plot the surface

% Note because of the way the base is defined, this is into the page
dCvec_prev = [0;0;1];
idCvec_prev = [1;0;0];

idC = -1i*dC;

% Normal vectors
v1 = [real(idC(1,:)); imag(idC(1,:)); reim(idC(2,:))];
v2 = [real(dC(1,:)); imag(dC(1,:)); reim(dC(2,:))];
nv = cross(v1, v2);
for k = 1:K
    nv(:,k) = nv(:,k)/norm(nv(:,k));
end

norm_prev = [0; 1; 0];

% figure(10)
plot3(real(C(1,:)), imag(C(1,:)), reim(C(2,:)), 'k'); hold all
    plot3(base(1,:), base(2,:), base(3,:), 'o-'); 
    
% We run through each step in the contour
for k = 1:K
    
    dCvec = [real(dC(1,k)); imag(dC(1,k)); reim(dC(2,k))];
    dCvec = dCvec/norm(dCvec);
    base = rotateObjectTwoVectors(base, dCvec_prev, dCvec);    
    
    idCvec_prev = rotateObjectTwoVectors(idCvec_prev, dCvec_prev, dCvec); 
    norm_prev = rotateObjectTwoVectors(norm_prev, dCvec_prev, dCvec); 
    dCvec_prev = dCvec;
    
    base = rotateObjectTwoVectors(base, norm_prev, nv(:,k));
    dCvec_prev = rotateObjectTwoVectors(dCvec_prev, norm_prev, nv(:,k));
    norm_prev = nv(:,k);
%     SUR(:,:,k) = base + repmat([real(C(1,k)); imag(C(1,k)); reim(C(2,k))], 1, npt);
    dot(norm_prev, dCvec_prev)
    
    
    
    % Rotate also along the surface
    idCvec = [real(idC(1,k)); imag(idC(1,k)); reim(idC(2,k))];
    idCvec = idCvec/norm(idCvec);
    
    %{
    plot3(base(1,:), base(2,:), base(3,:), 'o-');
    tmp1 = [[0;0;0], idCvec_prev];
    tmp2 = [[0;0;0], idCvec];
    plot3(tmp1(1,:), tmp1(2,:), tmp1(3,:))
    plot3(tmp2(1,:), tmp2(2,:), tmp2(3,:), 'w', 'LineWidth', 2)
    
    
    
    base = rotateObjectTwoVectors(base, idCvec_prev, idCvec);
    dCvec_prev = rotateObjectTwoVectors(dCvec_prev, idCvec_prev, idCvec);
    
    plot3(base(1,:), base(2,:), base(3,:), '*-');
    keyboard
    
    idCvec_prev = idCvec;
    dot(dCvec_prev, idCvec_prev)
    
    %}
    

    
    
    SUR(:,:,k) = base + repmat([real(C(1,k)); imag(C(1,k)); reim(C(2,k))], 1, npt);
%     plot3(SUR(1,:,k), SUR(2,:,k), SUR(3,:,k), 'o-', 'LineWidth', 2);
    pt = [real(C(1,k)); imag(C(1,k)); reim(C(2,k))];
%     tmn = [pt, pt + nv(:,k)];
    tmn = [pt, pt + norm_prev];
%     plot3(tmn(1,:), tmn(2,:), tmn(3,:), 'r', 'LineWidth', 1);
    
    test_tang = [pt, pt + dCvec];
%     plot3(test_tang(1,:), test_tang(2,:), test_tang(3,:), 'g', 'LineWidth', 1);
    
%     keyboard
    
end


X = squeeze(SUR(1,:,:));
Y = squeeze(SUR(2,:,:));
Z = squeeze(SUR(3,:,:));

% Fix any NaN's (This should not happen...)
xisnan = find(isnan(X(1,:)), 1);
xnotnan = find(~isnan(X(1,:)));
if ~isempty(xisnan)
    warning('NaN''s found');
    X = X(:,xnotnan);
    Y = Y(:,xnotnan);
    Z = Z(:,xnotnan);
    C = C(:,xnotnan);
end


surf(X,Y,Z)
axis equal
dragzoom;



