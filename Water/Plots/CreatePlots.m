%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Hamed Rahimian 
%           The Ohio State University 
%           January 01, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variation DIstance Constrained Water problem


n=200;
k=21;

scatter3(Frequency(:,2), Frequency(:,3), 100*Frequency(:,4)/21, 1/4*100*Frequency(:,4)/21, 'filled');
view(-30, 10);
xlabel('Avg. Supply');
ylabel('Avg. Demand');
zlabel('Frequency (%)');

