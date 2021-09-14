%% Calculate foot progression angle
function FPA = find_FPA(markers,whichLeg)
mkr1 = markers.mt2;
mkr2 = markers.cal;

footVec = mkr1 - mkr2;


for i = 1:length(footVec)
    % Horizontal plane
    FPA(i) = atand(footVec(2,i)/footVec(1,i));
end

% Decreasing FPA should cause KAM to decrease
% Flip sign for right leg
if strcmp(whichLeg,'RT'), FPA = -FPA; end

end