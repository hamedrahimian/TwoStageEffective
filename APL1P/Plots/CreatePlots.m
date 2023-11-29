%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Hamed Rahimian 
%           The Ohio State University 
%           January 01, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variation Distance Constrained APL1P 

noScen=[4, 5, 4, 4, 4];
n=1280;
A=zeros(n, 3);
for i=1:n;
    A(i,3)=100*Frequency(i,7)/Frequency(i,8);
    A(i,1)=(Frequency(i,2)+Frequency(i,3))/2;
    A(i,2)=(Frequency(i,4)+Frequency(i,5)+Frequency(i,6))/3;
end
scatter3(A(:,1), A(:,2), A(:,3), 1/4* A(:, 3), 'filled');
xlabel('Avg. Generator Availability');
ylabel('Avg. Demand');
zlabel('Frequency (%)');
view(-30, 10);




