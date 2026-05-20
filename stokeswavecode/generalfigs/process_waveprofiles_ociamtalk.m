% Processes wave profile data
%
% Nov 2015
% A version for the OCIAM talk

clear
close all

% Regular series
name = {...
    'series_N100_ep0.1.mat', ...
    'series_N100_ep0.22361.mat', ...
    'series_N100_ep0.31623.mat', ...
    'series_N100_ep0.44721.mat', ...
    'series_N100_ep0.54772.mat', ...
    'series_N100_ep0.63246.mat', ...
    'series_N100_ep0.70711.mat', ...
    'series_N100_ep0.7746.mat', ...
    'series_N100_ep0.83666.mat', ...
    'series_N100_ep0.89443.mat'};

p=0.1;
q = 0.4;
path = [p, p+q*1i, -p+q*1i, -p, p, p+q*1i, -p+q*1i, -p, p];

figure(1); clf(1); hold all;
for j = 1:length(name)
    load(['../../stokeswavedata/pubdata/profiles/', name{j}]);
    disp(['ep = ', num2str(solseries.ep)]);
    disp(['mu = ', num2str(solseries.mu, 12)]);
    
    [phi, z, zp] = get_profile_regularsol(solseries);
    
    % Phil modification to plot two periods
    x_add = [real(z) - 1, real(z)];
    y_add = [imag(z), imag(z)];
    
    plot(x_add, y_add - imag(z(1)), 'LineWidth', 2); 
    
    
    [f,g,gp,~,singloc] = gbern_contour(path,solseries);
    
    path_sing = [0, 1i*imag(singloc)-1i*1e-5];
%     path_sing = [0, singloc];
    [f2,g2,gp2,~,~] = gbern_contour(path_sing,solseries);

    % Then we calculate the y value at the singularity. g(f) = z(f) + z(-f)
    [z_sing_minus,~] = zfree(-singloc, solseries);
    z_sing(j) = g2(end) - z_sing_minus;

    f_sing(j) = singloc;
    H = imag(-zfree(0.5, solseries) + zfree(0, solseries));
    
    disp(['f_sing = ', num2str(f_sing(j), 16)]);
    disp(['z_sing = ', num2str(z_sing(j), 16)]);
    disp(['y_sing - y(-0.5)= ', num2str(imag(z_sing(j) - z(1)), 16)]);
    disp(['H = ', num2str(H, 16)]);
    disp('-----------------------------');
    
%     plot(0, imag(z_sing(j)) - imag(z(1)), 'ko');
    
    
end

% Steep series
name = {...
%     'series_steep_N120_ep0.94868.mat', ...
%     'series_steep_N120_ep0.96954.mat', ...
    'series_steep_N120_ep0.99499.mat'};

for j = 1:length(name)
    load(['../../stokeswavedata/pubdata/profiles/', name{j}]);
    disp(['ep = ', num2str(steepsolseries.ep)]);
    disp(['mu = ', num2str(steepsolseries.mu, 6)]);
    [phi, z, zp] = get_profile_steepsol(steepsolseries);
    
    x_add = [real(z) - 1, real(z)];
    y_add = [imag(z), imag(z)];
    plot(x_add, y_add - imag(z(1)), 'LineWidth', 2);
        
end

xlabel('x');
ylabel('y - y(-0.5)');
