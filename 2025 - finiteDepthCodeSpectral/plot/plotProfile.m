close all


data = load("energy=0.1.mat")
obj = figure(1);
theme(obj,'light');
hold on
plot3(data.xi,data.AmpV.*ones(1024,5),data.Yout,Color='k')
plot3(data.xi,data.AmpV.*ones(1024,5),zeros(1024,5),Color='k')


data = load("energy=0.2.mat")
obj = figure(2);
theme(obj,'light');
hold on
plot3(data.xi,data.AmpV.*ones(1024,5),data.Yout,Color='k')
plot3(data.xi,data.AmpV.*ones(1024,5),zeros(1024,5),Color='k')