function [P, N, pind, nind]= perc_unperc_overlap_check_function(P, N, pind, nind)
% This function takes the percolation matrix and unpercolated matrix and
% checks for overlap

% check for overlaping
[skip_int, intp, intn] = intersect(P(1:pind-1,:),N(1:nind-1,1:3),'rows');
overlap_vox_id = unique(N(intn, 4));
% if vesicle overlaps
if ~isempty(overlap_vox_id)
    for ii = 1:length(overlap_vox_id)
        if overlap_vox_id(ii) == 0
            continue
        end
        new_perc_ind = find(N(:,4) == overlap_vox_id(ii));
        lnp = length(new_perc_ind);
        % Add to perc matrix
        P(pind:pind+lnp-1,:) = N(new_perc_ind, 1:3);
        pind = pind+lnp;
        % Remove from nonperc matrix
        N(new_perc_ind, :) = [];
        nind = nind -lnp;
    end
end