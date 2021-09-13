function store_wrt_pelvis = wrtPelvis(input, pelvis)
store_wrt_pelvis = cell(12,10);

for i = 1:1:12
    for j = 1:1:10
        if i == 3||i==4||i==6||i==12 % if R leg, negative in ML direction
            store_wrt_pelvis{i,j} = [input{i,j}(1,:) - pelvis{i,j}(1,:); -input{i,j}(2,:) + pelvis{i,j}(2,:); input{i,j}(3,:) - pelvis{i,j}(3,:)];
            %store_wrt_pelvis{i,j} = [input{i,j}(1,:) - pelvis{i,j}(1,:); input{i,j}(2,:) - pelvis{i,j}(2,:); input{i,j}(3,:) - pelvis{i,j}(3,:)];
        
        else
           store_wrt_pelvis{i,j} = [input{i,j}(1,:) - pelvis{i,j}(1,:); input{i,j}(2,:) - pelvis{i,j}(2,:); input{i,j}(3,:) - pelvis{i,j}(3,:)];
        end
    end
end

end