figure(1)
contourf(XI,ETA,real(YmatFull),50)

figure(2)
plot3(XI',ETA',real(YmatFull)',Color='Blue');
hold on
plot3(xi_mesh,xi_mesh*0,real(Y0),Color='red')
grid off
hold off

figure(3)
p = plot3(XI',ETA',real(YmatFull)',Color='black');
hold on
plot3(xi_mesh,xi_mesh*0,real(Y0),Color='red')
grid off
view(0,0)
xlabel('\xi')
ylabel('real(Y)')
%ylim([0.8,0.15])
hold off