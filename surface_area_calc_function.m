function [san, sap]= surface_area_calc_function(P, N, Vsa, pind, nind, vind)
% This function calculates the number of interfaces of vesicles that are
% touching matrix

san = 0;
sap = 0;

% Determine which surface area voxesl intersect with perc matrix
[int] = intersect(P(1:nind-1,:), Vsa(1:vind-1,:),'rows');

% loop through percolation matrix
if pind>100
    % create surface matrix if all voxels surround all surface voxels
    surface_mat = [[int(:,1)+1, int(:,2:3)]; ...
        [int(:,1)-1, int(:,2:3)]; ...
        [int(:,1), int(:,2)+1, int(:,3)]; ...
        [int(:,1), int(:,2)-1, int(:,3)]; ...
        [int(:,1:2), int(:,3)+1]; ...
        [int(:,1:2), int(:,3)-1]];
    % remove overlaping voxels
    surface_mat = unique(surface_mat, 'rows');
    [int_suf] = intersect(P(1:pind,:), surface_mat,'rows');
    
    sap = length(surface_mat) - length(int_suf);
end

% now check unpercolated matrix
% Determine which surface area voxesl intersect with nonperc matrix
[int] = intersect(N(1:nind-1,1:3),Vsa(1:vind-1,:),'rows');

if nind > 100
    % create surface matrix if all voxels surround all surface voxels
    surface_mat = [[int(:,1)+1, int(:,2:3)]; ...
        [int(:,1)-1, int(:,2:3)]; ...
        [int(:,1), int(:,2)+1, int(:,3)]; ...
        [int(:,1), int(:,2)-1, int(:,3)]; ...
        [int(:,1:2), int(:,3)+1]; ...
        [int(:,1:2), int(:,3)-1]];
    
    % remove overlaping voxels
    surface_mat = unique(surface_mat,  'rows');
    [int_suf] = intersect(N(1:nind,1:3), surface_mat,'rows');
    
    san = length(surface_mat) - length(int_suf);
end

%%
%
% % tic
% % Determine which surface area voxesl intersect with perc matrix
% [int, intp, intvp] = intersect(P(1:nind,:), Vsa(1:vind,:),'rows');
% rem_intvp = zeros(length(intvp),1);
% % loop through percolation matrix
% if pind>100
%     for i=2:length(int)
%         current_vox = int(i,:);
%         surface_mat = [current_vox(1)+1, current_vox(2:3); ...
%             current_vox(1)-1, current_vox(2:3); ...
%             current_vox(1), current_vox(2)+1, current_vox(3); ...
%             current_vox(1), current_vox(2)-1, current_vox(3); ...
%             current_vox(1:2), current_vox(3)+1; ...
%             current_vox(1:2), current_vox(3)-1];
%
%         [~, ~, c] = intersect(P(i:pind,:),surface_mat,'rows');
%         number_faces = length(c);
%         if number_faces > 0
%             rem_intvp(i) = 1;
%             sap = sap + number_faces;
%         end
%     end
% end
%
% % remove Vsa values that don't touch matrix
% Vsa(intvp(rem_intvp==0),:) = [];
% vind = vind - length(find(rem_intvp==0));
% %% loop through unpercolated matrix
% % Determine which surface area voxesl intersect with nonperc matrix
% [int, intn, intvn] = intersect(N(1:nind,1:3),Vsa(1:vind,:),'rows');
% rem_intvn = zeros(length(intvn),1);
% if nind > 100
%     for i=2:length(int)
%         current_vox = int(i,:);
%
%         surface_mat = [current_vox(1)+1, current_vox(2:3); ...
%             current_vox(1)-1, current_vox(2:3); ...
%             current_vox(1), current_vox(2)+1, current_vox(3); ...
%             current_vox(1), current_vox(2)-1, current_vox(3); ...
%             current_vox(1:2), current_vox(3)+1; ...
%             current_vox(1:2), current_vox(3)-1];
%
%         [~, ~, c] = intersect(N(i:nind,1:3),surface_mat,'rows');
%
%         number_faces = length(c);
%
%         if number_faces > 0
%             rem_intvn(i) = 1;
%             san = san  + number_faces;
%         end
%     end
% end
% % remove Vsa values that don't touch matrix
% Vsa(intvn(rem_intvn==0),:) = [];
% vind = vind - length(find(rem_intvn==0));
% % toc