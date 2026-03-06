%% =========================================================
%  HELPER FUNCTION
%% =========================================================

function seg = kpath_segment(k_start, k_end, n_pts)
% Returns n_pts row vectors linearly interpolated from k_start to k_end
    t   = linspace(0, 1, n_pts)';          % (n_pts x 1)
    dk  = k_end(:)' - k_start(:)';         % ensure row (1 x 2)
    seg = k_start(:)' + t * dk;            % (n_pts x 2)
end