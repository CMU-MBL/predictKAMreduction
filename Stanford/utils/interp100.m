function resized = interp100(input_store)

% initialize output
resized = cell(12,10);

for i = 1:1:12
    for j = 1:1:10
        % break down into x,y,z
        inp_x = input_store{i,j}(1,:);
        inp_y = input_store{i,j}(2,:);
        inp_z = input_store{i,j}(3,:);
        
        % resample to 100
        x = interp1(linspace(1, length(inp_x), length(inp_x)), inp_x, linspace(1, length(inp_x), 100));
        y = interp1(linspace(1, length(inp_y), length(inp_y)), inp_y, linspace(1, length(inp_y), 100));
        z = interp1(linspace(1, length(inp_z), length(inp_z)), inp_z, linspace(1, length(inp_z), 100));
        
        % save to resized cell
        resized{i,j} = [x; y; z];
    end
end

end