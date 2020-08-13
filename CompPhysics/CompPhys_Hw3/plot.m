function plot
%%PROBLEM 7.2
load('sunspots.mat')
plot(sunspots(:,1),sunspots(:,2),'k-');
xlabel('Time (months)');
ylabel('Sunspot Number');
%%
end

