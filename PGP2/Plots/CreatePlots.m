%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Hamed Rahimian 
%           The Ohio State University 
%           January 01, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variation Distance Constrained PGP2 


noScen=[9, 8, 8];
n=576;

A=zeros(n, 3);
for i=1:n;
    A(i,3)=100*Frequency(i,5)/Frequency(i,6);
    A(i,1)=Frequency(i,2);
    A(i,2)=(Frequency(i,3)+Frequency(i,4))/2;
end
scatter3(A(:,1), A(:,2), A(:,3), 1/4*A(:,3), 'filled');
view(20, 10);
grid on;
xlabel('Demand 1');
ylabel('Avg. Demand 2 and 3');
zlabel('Frequency (%)');
figure;

B=zeros(n, 3);
for i=1:n;
    B(i,3)=100*Frequency(i,5)/Frequency(i,6);
    B(i,1)=Frequency(i,3);
    B(i,2)=(Frequency(i,2)+Frequency(i,4))/2;
end
scatter3(B(:,1), B(:,2), B(:,3), 1/4*B(:,3), 'filled');
view(20, 10);
grid on;
xlabel('Demand 2');
ylabel('Avg. Demand 1 and 3');
zlabel('Frequency (%)');
figure;
C=zeros(n, 3);
for i=1:n;
    C(i,3)=100*Frequency(i,5)/Frequency(i,6);
    C(i,1)=Frequency(i,4);
    C(i,2)=(Frequency(i,2)+Frequency(i,3))/2;
end
scatter3(C(:,1), C(:,2), C(:,3), 1/4*C(:,3), 'filled');
view(20, 10);
grid on;
xlabel('Demand 3');
ylabel('Avg. Demand 1 and 2');
zlabel('Frequency (%)');

