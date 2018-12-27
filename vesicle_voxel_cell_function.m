function [C, lc, Csa, lv]= vesicle_voxel_cell_function(rad_vox, nv, vol_dim)
% This function takes the vesicle center coordinate and radius and
% determines voxel index (x,y,z) of all voxels in vesicle

% find bounds of box surrounding voxel
lowx = round(nv(1)-rad_vox);
if lowx < 1
    lowx = 1;
end
lowy = round(nv(2)-rad_vox);
if lowy < 1
    lowy = 1;
end
lowz = round(nv(3)-rad_vox);
if lowz < 1
    lowz = 1;
end

hx = round(nv(1)+rad_vox);
if hx > vol_dim
    hx = vol_dim;
end
hy = round(nv(2)+rad_vox);
if hy > vol_dim
    hy = vol_dim;
end
hz = round(nv(3)+rad_vox);
if hz > vol_dim
    hz = vol_dim;
end

% initialize current vesicle coordinate matrix
C = uint16(zeros(round(rad_vox*5/3)^3,3));
% initialize current vesicle surface area coordinate matrix
Csa = uint16(zeros(round(rad_vox*5/3)^3,3));
% Search vesicles in box to find which ones are in sphere
m = 1;
mm = 1;
for i = lowx:hx
    for j = lowy:hy
        for k = lowz:hz
            % distance of current voxel from vesicle center
            vox_dist = sqrt((i-nv(1))^2 + (j-nv(2))^2 + (k-nv(3))^2);
            % if voxel is within the radius of the vesicle
            if vox_dist < rad_vox
                C(m,:) = [i,j,k];
                m = m+1;
                % if the voxel is near edge of vesicle then save it
                if vox_dist > rad_vox-1.5
                    Csa(mm,:) = [i,j,k];
                    mm = mm+1;
                end
            end
        end
    end
end

% remove unneeded preallocated rows
C = unique(C, 'rows');
C= C(2:end,:);
% length C
lc = length(C(:,1));
% remove unneeded preallocated rows
Csa = unique(Csa, 'rows');
Csa= Csa(2:end,:);
% length C
lv = length(Csa(:,1));

% meanj = mean([lowy,hy]);
%  ind_vox = find(Csa(:,2)== round(45));
%  center_slice = zeros(hx-lowx, hz-lowz);
%         P_image = Csa(ind_vox,[1,3]);
%         for p=1:length(P_image)
%             center_slice(P_image(p,1),P_image(p,2)) = 2;
%         end
%         
%         figure
%         h = imagesc(center_slice);
% %         set(h,'alphadata',center_slice ~=0)
%         axis equal
%         axis tight