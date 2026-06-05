%% A script to plot real solution and location of primary singularity


obj = figure(6);
theme(obj,'light');
hold on
plot(E,imag(Svec(end)),'o',Color='black')
plot(xi,1.4*(Yreal-min(Yreal)),Color='black')
text(0.3,imag(Svec(end)),['E=', num2str(E, '%.1e')])
plot([0,0.3],[imag(Svec(end)),imag(Svec(end))],'--',Color='black')

text(-0.3,1.4*(max(Yreal)-min(Yreal)),['E=', num2str(E, '%.1e')],'HorizontalAlignment','right')
plot([0,-0.3],1.4*[max(Yreal)-min(Yreal),max(Yreal)-min(Yreal)],'--',Color='black')
hold off
ylabel('Im(f)')
xlabel('X')
title(['M=', num2str(MV)])