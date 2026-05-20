clear 
close all

% Singularity location data for ep=0.1, ep=0.2. Ep is sufficiently small
% to use the boundary integral for singularities (no need for steep series)
% Using singularity notation from explanation.jpg. E.g. B on BDE means the
% B-type singularity on the sheet BDE.

eps = [0.1,0.15,0.2];

for k=1:length(eps)

% Load in the solution
ep = eps(k);
loadpath = ['../../stokeswavedata/series_N100_ep', num2str(ep), '.mat'];
F = load(loadpath);
solseries = F.solseries;

%{
% a on 0 
path = get_path('a-0'); 
sing_start=1; 
[~, singloc(k,1)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start);
singloc(k,1) = -singloc(k,1); % Because of path orientation

% A on a
path = get_path('A-a');
sing_start=6;
[~, singloc(k,2)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start);
singloc(k,2) = -singloc(k,2); 

% B^ on aA
path = get_path('B^-aA');
sing_start=8; 
[~, singloc(k,3)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start);
singloc(k,3) = singloc(k,3); 

% B on aA
path = get_path('B-aA'); 
sing_start=9;
[~, singloc(k,4)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start);
singloc(k,4) = -singloc(k,4); 

% C on aAaA
path = get_path('C-aAaA');
sing_start=18; 
[~, singloc(k,5)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start);
singloc(k,5) = singloc(k,5); 

% C^ on aAaA
path = get_path('C^-aAaA');
sing_start=18; 
[~, singloc(k,6)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start);
singloc(k,6) = singloc(k,6); 

% a on aA
path = get_path('a-aA');
sing_start=11; 
[~, singloc(k,7)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start);
singloc(k,7) = singloc(k,7); 

% a on aAB^
path = get_path('a-aAB^');
sing_start=12; 
[~, singloc(k,8)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start);
singloc(k,8) = singloc(k,8); 
%}

% b on aAB^a
path = get_path('b-aAB^a');
sing_start=14; 
[~, singloc(k,9)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start);
singloc(k,9) = -singloc(k,9); 

%{
% b^ on aAB^a
path = get_path('b^-aAB^a');
sing_start=16; 
[~, singloc(k,10)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start);
singloc(k,10) = singloc(k,10); 
%}

end

% Plot stuff
plot(real(singloc(1,:)), imag(singloc(1,:)), 'b*');
hold on
plot(real(singloc(2,:)), imag(singloc(2,:)), 'r*');
plot(real(singloc(3,:)), imag(singloc(3,:)), 'g*');
%plot(real(singloc(4,:)), imag(singloc(4,:)), 'y*');
%plot(real(singloc(5,:)), imag(singloc(5,:)), 'c*');


%save('../singularitydat/singularitydat.mat', 'singloc');