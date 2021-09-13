function [pelvis_center_baseline, pelvis_center_trial] = find_PelvisCenter(i, bmarkers, tmarkers)


    % load ASIS and PSIS vals for baseline and trial
    %           LASIS baseline
    lasi_x_b =  bmarkers.lasi(1, :);
    lasi_y_b = bmarkers.lasi(2, :);
    lasi_z_b = bmarkers.lasi(3, :);
    
    %           RASIS baseline
    rasi_x_b =  bmarkers.rasi(1, :);
    rasi_y_b = bmarkers.rasi(2, :);
    rasi_z_b = bmarkers.rasi(3, :);
    
    %           LPSIS baseline
    lpsi_x_b =  bmarkers.lpsi(1, :);
    lpsi_y_b = bmarkers.lpsi(2, :);
    lpsi_z_b = bmarkers.lpsi(3, :);
    
    %           RPSIS baseline
    rpsi_x_b =  bmarkers.rpsi(1, :);
    rpsi_y_b = bmarkers.rpsi(2, :);
    rpsi_z_b = bmarkers.rpsi(3, :);
    
    %           LASIS trial
    lasi_x_t = tmarkers.lasi(1,:);
    lasi_y_t = tmarkers.lasi(2,:);
    lasi_z_t = tmarkers.lasi(3,:);
    
    %           RASIS trial
    rasi_x_t = tmarkers.rasi(1,:);
    rasi_y_t = tmarkers.rasi(2,:);
    rasi_z_t = tmarkers.rasi(3,:);
    
    %           LPSIS trial
    lpsi_x_t = tmarkers.lpsi(1,:);
    lpsi_y_t = tmarkers.lpsi(2,:);
    lpsi_z_t = tmarkers.lpsi(3,:);
    
    %           RPSIS trial
    rpsi_x_t = tmarkers.rpsi(1,:);
    rpsi_y_t = tmarkers.rpsi(2,:);
    rpsi_z_t = tmarkers.rpsi(3,:);
    
    % s6 LASIS trial: remove marker artefact noise at end
    if i == 6
        lasi_x_t =  tmarkers.lasi(1, 1:length(lpsi_x_t));
        lasi_y_t = tmarkers.lasi(2, 1:length(lpsi_y_t));
        lasi_z_t = tmarkers.lasi(3, 1:length(lpsi_z_t));
    end
    
   % s11 RASIS baseline and trial: remove marker artefact noise at end
    if i == 11
        rasi_x_b =  bmarkers.rasi(1, 1:length(lasi_x_b));
        rasi_y_b =  bmarkers.rasi(2, 1:length(lasi_y_b));
        rasi_z_b =  bmarkers.rasi(3, 1:length(lasi_z_b));
        
        rasi_x_t = tmarkers.rasi(1, 1:length(lasi_x_t));
        rasi_y_t = tmarkers.rasi(2, 1:length(lasi_y_t));
        rasi_z_t = tmarkers.rasi(3, 1:length(lasi_z_t));
    end
    
    % calculate centroids of 4 subtriangles
    %       T1: LASI, RASI, LPSI
    %           baseline
    t1_x_b = (lasi_x_b + rasi_x_b + lpsi_x_b)/3;
    t1_y_b = (lasi_y_b + rasi_y_b + lpsi_y_b)/3;
    t1_z_b = (lasi_z_b + rasi_z_b + lpsi_z_b)/3;
    
    %           trial
    t1_x_t = (lasi_x_t + rasi_x_t + lpsi_x_t)/3;
    t1_y_t = (lasi_y_t + rasi_y_t + lpsi_y_t)/3;
    t1_z_t = (lasi_z_t + rasi_z_t + lpsi_z_t)/3;
    
    %       T2: LASI, RASI, RPSI
    %           baseline
    t2_x_b = (lasi_x_b + rasi_x_b + rpsi_x_b)/3;
    t2_y_b = (lasi_y_b + rasi_y_b + rpsi_y_b)/3;
    t2_z_b = (lasi_z_b + rasi_z_b + rpsi_z_b)/3;
    
    %           trial
    t2_x_t = (lasi_x_t + rasi_x_t + rpsi_x_t)/3;
    t2_y_t = (lasi_y_t + rasi_y_t + rpsi_y_t)/3;
    t2_z_t = (lasi_z_t + rasi_z_t + rpsi_z_t)/3;
    
    %       T3: RASI, RPSI, LPSI
    %           baseline
    t3_x_b = (lpsi_x_b + rasi_x_b + rpsi_x_b)/3;
    t3_y_b = (lpsi_y_b + rasi_y_b + rpsi_y_b)/3;
    t3_z_b = (lpsi_z_b + rasi_z_b + rpsi_z_b)/3;
    
    %           trial
    t3_x_t = (lpsi_x_t + rasi_x_t + rpsi_x_t)/3;
    t3_y_t = (lpsi_y_t + rasi_y_t + rpsi_y_t)/3;
    t3_z_t = (lpsi_z_t + rasi_z_t + rpsi_z_t)/3; 
    
    %       T4: RPSI, LPSI, LASI
    %           baseline
    t4_x_b = (lpsi_x_b + lasi_x_b + rpsi_x_b)/3;
    t4_y_b = (lpsi_y_b + lasi_y_b + rpsi_y_b)/3;
    t4_z_b = (lpsi_z_b + lasi_z_b + rpsi_z_b)/3;
    
    %           trial
    t4_x_t = (lpsi_x_t + lasi_x_t + rpsi_x_t)/3;
    t4_y_t = (lpsi_y_t + lasi_y_t + rpsi_y_t)/3;
    t4_z_t = (lpsi_z_t + lasi_z_t + rpsi_z_t)/3; 
    
    
    % (system of equations approach wouldn't work because the 2 lines don't necessarily intersect in 3D space)
    % instead: take average of 4 centroids to find x, y, z coordinates
    
    % baseline
    x_b = (t1_x_b + t2_x_b + t3_x_b + t4_x_b)/4;
    y_b = (t1_y_b + t2_y_b + t3_y_b + t4_y_b)/4;
    z_b = (t1_z_b + t2_z_b + t3_z_b + t4_z_b)/4;
    
    % trial
    x_t = (t1_x_t + t2_x_t + t3_x_t + t4_x_t)/4;
    y_t = (t1_y_t + t2_y_t + t3_y_t + t4_y_t)/4;
    z_t = (t1_z_t + t2_z_t + t3_z_t + t4_z_t)/4;
    
    % OUTPUT: 
    % (x, y, z) form:
    pelvis_center_baseline{i} = [x_b; y_b; z_b];
    pelvis_center_trial{i} = [x_t; y_t; z_t];

end