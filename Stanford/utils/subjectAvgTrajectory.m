function [avg_x, avg_y, avg_z] = subjectAvgTrajectory (input_store)
% initialize 
avg_x = zeros(12,100);
avg_y = zeros(12,100);
avg_z = zeros(12,100);

for ii = 1:1:12 % loop through every subject
    sumx = 0;
    sumy = 0;
    sumz = 0;
    for j = 1:1:10 % loop through every step and compute running sum
        sumx = sumx + input_store{ii,j}(1,:);
        sumy = sumy + input_store{ii,j}(2,:);
        sumz = sumz + input_store{ii,j}(3,:);
    end
    avg_x(ii,:) = sumx/10; % save parameter average for subject #ii
    avg_y(ii,:) = sumy/10;
    avg_z(ii,:) = sumz/10;

end