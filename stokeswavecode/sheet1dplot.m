close all
clear all

F = load('data/recdat_rec9_ep0.8.mat');
recdat = F.recdat;

recnum = length(recdat);

%{
x=[];
y=[];
for k=1:recnum
    x = [x; transpose(linspace(k-1,k,length(recdat(k).z)))];
    y = [y; recdat(k).z];
end
%}

%recnum should be odd
x=[]; y=[];

for k=1:recnum
    x = [x; k];
    y = [y; recdat(k).z(end)];
end

figure(5);
hold on
plot(x,real(y),'o');
plot(x,imag(y),'o');

