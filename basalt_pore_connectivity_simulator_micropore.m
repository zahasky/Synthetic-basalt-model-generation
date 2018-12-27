% basalt_pore_connectivity_simulator
% Christopher Zahasky
% 12/15/2017 - updated for web 12/27/2018
% This is an adapted version of 'basalt_pore_connectivity_simulator.m' 
% for the generating synthetic basalt models with microporosity. 
% This works by iteratively added vesicles randomly to a 3D
% volume and calculting effective porosity as vesicles are added. NOTE THAT
% SURFACE AREA CALCULATIONS ARE UNTESTED as there was not enough
% information about microporous surface area to properly benchmark these
% types of simulation results.

clear all
close all

% simluator porosity at which simulation stops
phi_stop_threshold = 0.50;
% voxel coordinate spacing (for voxel_coord_matrix_function)
coord_spacing = 3;
% grid cell side length, add multiple values to run batch of simulations
% with different model size
VOLD = [120];
micro_pore_thick = 2;

% loop for running different size models
for v = 1:length(VOLD)
    v
    vol_dim = VOLD(v); 
    total_vol = vol_dim^3;
    % Calculate the surface area every few added voxels (to save time)
    sa_calc_interval = total_vol*0.001;
    min_vesicle_rad_vox =  2;
    max_vesicle_rad_vox = 9;
    
    % loop for running specified number of simulations with give model size
    number_of_simulations = 50;
    for r = 1:number_of_simulations
        % output simulation number to track progress
        r
        tic
        % Set voxel coordinates
        [seed_coord]= voxel_coord_matrix_function(vol_dim, coord_spacing);
        
        % itit size of variables
        init_size = vol_dim;
        % initialize variables
        phi = zeros(init_size,1);
        phi_eff = zeros(init_size,1);
        phi_path = zeros(init_size,1);
        sa_perc = zeros(init_size,1);
        sa_non = zeros(init_size,1);
        
        % preallocate 'percolated' coordinate matrix [x y z]
        pind = 1;
        pind_last = 1;
        P = uint16(zeros(vol_dim^3, 3));
        % preallocate 'percolated' micropore coordinate matrix [x y z]
        pindmp = 1;
%         pind_lastmp = 1;
        Pmp = uint16(zeros(vol_dim^3, 3));
        % preallocate 'unpercolated' matrix
        nind = 1;
        nind_last = 1;
        N = uint16(zeros(vol_dim^3, 4));
        nn =uint16(1);
        n=uint16(1);
        
        % generate random radius matrix type 1 = rand, 2 = normal, 3 = log, 4 =
        % bimodal
        mat_length = 5000;
        [rad_vox_mat] = vesicle_rad_function(2, ...
            min_vesicle_rad_vox, max_vesicle_rad_vox, mat_length);

        % vesicle surface area matrix
        Vsa = uint16(zeros(vol_dim^3, 3));
        vind = 1;
        jet_range = 100;
        custom_jet = [0 0 0; jet(jet_range)];
        
        while phi < phi_stop_threshold
            % In order to speed up algorithm and prevent repetative seeding in same
            % location use randomly draw number from seed_coord matrix
            % new vesicle center coord
            seed_ind = round(rand(1)*(length(seed_coord)-1) + 1);
            nv = seed_coord(seed_ind,:);
            seed_coord(seed_ind,:) = [];
            % Select random radius
            rad_vox = rad_vox_mat(rem(n, mat_length)+1);
            % Find all voxels in vesicle
            [C, lc, Csa, lv] = vesicle_voxel_cell_function(rad_vox, nv, vol_dim);
            % Find all voxels in microporous vesicle
            [Cmp, lcmp] = vesicle_voxel_cell_function(rad_vox+micro_pore_thick, nv, vol_dim);
            % update surface area matrix
            Vsa(vind:vind+lv-1,:) = Csa;
            vind = vind + lv;
            
            % Now determine if vesicle crosses inlet
            inlet_intersect = find(Cmp(:,3)==1);
            
            % Look for any overlap between new vesicle voxels and percolated
            % voxels
            [int, intp, intc] = intersect(Pmp(1:pindmp,:),Cmp,'rows');
            % Look for any overlap between new vesicle voxels and percolated
            % voxels
            [non_int, intn, intcc] = intersect(N(1:nind,1:3),Cmp,'rows');
            % if voxels overlap save as percolated
            if ~isempty(int)
%                 C(intc,:) = [];
%                 lc = length(C(:,1));
%                 P(pind:pind+lc-1,:) = C;
%                 pind = pind+lc;
                
                Cmp(intc,:) = [];
                lcmp = length(Cmp(:,1));
                Pmp(pindmp:pindmp+lcmp-1,:) = Cmp;
                pindmp = pindmp+lcmp;
                
                % If there is overlap then select all nonpercolated voxels
                % with same voxel ID and convert to percolated matrix
                if ~isempty(non_int)
                    overlap_vox_id = N(intn, 4);
                    overlap_vox_id = unique(overlap_vox_id);
                    % find all voxels with overlap_vox_id
                    for pp=1:length(overlap_vox_id)
                        new_perc_ind = find(N(1:nind,4) == overlap_vox_id(pp));
                        lnp = length(new_perc_ind);
                        % Add to perc matrix
                        Pmp(pindmp:pindmp+lnp-1,:) = N(new_perc_ind, 1:3);
                        pindmp = pindmp+lnp;
                        % Remove from nonperc matrix
                        N(new_perc_ind, :) = [];
                        nind = nind -lnp;
                    end
                end
                
                % if no overlap with percolated voxels but there it is located at
                % inlet
            elseif ~isempty(inlet_intersect)
                Pmp(pindmp:pindmp+lcmp-1,:) = Cmp;
                pindmp = pindmp+lcmp;
                
                % if voxels don't overlap then save as unpercolated
            else
                % If there is overlap then set voxel id to be the same as
                if ~isempty(non_int)
                    overlap_vox_id = unique(N(intn, 4));
                    % if new vesicle connects two sets of unconnected vesicle
                    % assign same overlap id
                    if length(overlap_vox_id)>1
                        ddd = 1;
                    end
                    for pp=1:length(overlap_vox_id)
                        n_joined_ind = find(N(1:nind,4) == overlap_vox_id(pp));
                        N(n_joined_ind, 4) = overlap_vox_id(1);
                    end
                    vox_id_vec = uint16(ones(lcmp,1)).*overlap_vox_id(1);
                else
                    vox_id_vec = uint16(ones(lcmp,1)).*n;
                end
                % remove overlapping voxels
                Cmp(intcc,:) = [];
                lcmp = length(Cmp(:,1));
                N(nind:nind+lcmp-1,:) = [Cmp, vox_id_vec(1:lcmp)];
                nind = nind+lcmp;
            end
            
            %% Total porosity tracker
            % Look for any overlap between new vesicle voxels and previous
            % voxels (don't care about percolation)
            [int, intp, intc] = intersect(P(1:pind,:),C,'rows');
            % remove overlapped voxels
            C(intc,:) = [];
            lc = length(C(:,1));
            P(pind:pind+lc-1,:) = C;
            pind = pind+lc;
              
            
            % If ther has been some amount of difference in voxels added or lost
            % from matrix then calculate porosity and surface area
            if abs(nind-nind_last) > sa_calc_interval || abs(pind-pind_last) > sa_calc_interval
                nind_last = nind;
                pind_last = pindmp;
                % check for overlaping
                [Pmp, N, pindmp, nind]= perc_unperc_overlap_check_function(Pmp, N, pindmp, nind);
                % check for duplicates in vesicle matrix
                A = unique(P(1:pind),'rows');
                if length(A) ~= pind
                    warning('There are duplicates in P matrix')
                end
                % porosity calculated by vesicle volume without micropores
                phi(nn) = (pind)/total_vol;
%                 phi(nn) = (nind + pind)/total_vol;
                phi_eff(nn) = pindmp/total_vol;
                if phi(nn)<phi_eff(nn)
                    phi(nn) = phi_eff(nn);
                end
                if ~isempty(find(Pmp(1:pindmp,3)==vol_dim))
                    phi_path(nn) = phi_eff(nn);
                else
                    phi_path(nn) = 0;
                end
                
                [san, sap]= surface_area_calc_function(Pmp, N, Vsa, pindmp, nind, vind);
                sa_non(nn) = san;
                sa_perc(nn) = sap;
                nn = nn+1;
            end
            
            % plot a slice every 100 iterations
            if mod(n, 50)==0
                figure(1)
                % Find voxels that intersect center of volume
                ind_vox = find(P(1:pind,1)==round(vol_dim/2));
                P_image = P(ind_vox,[2:3]);
                center_slice = zeros(vol_dim, vol_dim);
                for p=1:length(P_image(:,1))
                    center_slice(P_image(p,1),P_image(p,2))=1;
                end
                
                ind_vox = find(N(1:nind,1)==round(vol_dim/2));
                P_image = N(ind_vox,[2:4]);
                for p=1:length(P_image(:,1))
                    center_slice(P_image(p,1),P_image(p,2)) = rem(P_image(p,3), jet_range)+1;
                end
                figure(1)
                h = imagesc(center_slice);
                set(h,'alphadata',center_slice ~=0)
                axis equal
                axis tight
                colormap(custom_jet)
                caxis([1 jet_range])
                title(['n index = ', num2str(n), ', phi total = ', ...
                    num2str(phi(nn-1)), ', phi effective = ', ...
                    num2str(phi_eff(nn-1))])
                drawnow
            end
            
            n = n+1;
        end
        sa_non(phi==0) = [];
        sa_perc(phi==0) = [];
        phi_eff(phi==0) = [];
        phi_path(phi==0) = [];
        phi(phi==0) = [];
        
        comp_data{r} = [phi, phi_eff, phi_path, sa_non, sa_perc];
        
        figure(3)
        hold on
        plot(phi, phi_eff, 'b')
        plot(phi, phi_path, 'r')
        plot(phi, phi, 'k')
        xlabel('Total Porosity')
        legend('Effective porosity', 'Connected porosity')
        drawnow
        toc
    end
    % save data
%     save(['normal_50sim_microt_2_vol_dim_', num2str(vol_dim)], ...
%         'comp_data', 'min_vesicle_rad_vox', 'max_vesicle_rad_vox', 'vol_dim', 'rad_vox_mat')
    clear comp_data
end
