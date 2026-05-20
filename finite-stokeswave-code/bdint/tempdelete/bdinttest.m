
clear
close all

N = 201; % Should be odd
ep = 0.4; % Amplitude
r0 = 0; % Depth parameter

name = ['./bdint_data/', 'bdint','_r', num2str(r0),'_ep', num2str(ep),'.mat'];
loaded_data = load(name);
keyboard
zp = loaded_data.zp;

delta_guess = real(zp);
beta_guess = imag(zp);

x = loaded_data.x;
y = loaded_data.y;
x2 = loaded_data.x2;
y2 = loaded_data.y2;


figure(2);
plot(x,y, '-'); hold all
plot(x2,y2, 'o');
axis normal

legend('Boundary integral', 'Series truncation');