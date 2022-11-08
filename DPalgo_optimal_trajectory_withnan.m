function projection_points_in_DSP = DPalgo_optimal_trajectory_withnan(DS_plane_feat, baby_data_feat_withnan, radius)
% To find the optimal projection of the data points onto the model
% [K^+]_{Buffer}-[O_2]_{Buffer} plane.
%
% Previous to this step, we perform feature normalization (for both DS_plane
% and baby_data) by subtracting each feature value from its global mean and
% dividing by the global standard deviation (z-scores).
%
% Input:
% DS_plane_feat : N x M x O matrix 
% baby_data_feat_withnan : P x O matrix
% radius : 1x1 scalar
%
% Output:
% projection_points_in_DSP : P x 2 vector
%


% filtering out rows with nan from baby data
baby_data = filtering_finite_rows(baby_data_feat_withnan);

N = size(DS_plane_feat,1); % length of [K^+]_{Buffer} dimension
M = size(DS_plane_feat,2); % length of [O_2]_{Buffer} dimension
O = size(DS_plane_feat,3); % num of features or burst-statistics

P = size(baby_data,1); % num of baby_points to be projected

OO = size(baby_data,2);

if O~=OO
    error('Error: The number DSP featurs and baby_data features should be same'); 
end


if P == 0 
    projection_points_in_DSP = nan;
    return
end


ppscore_single_row = zeros(1, N*M);

direction = zeros(P, N*M);


% [Filling the first row (N*M x 1)--- initial step]
% Filling the first row of ppscore matrix with the LN-norm of the first 
% data-point with all the points in the DS plane. This represents all the 
% possible starting points of the optimal projection. Similar to what we 
% had in the modified greedy approach.

for nmi = 1:N*M
    [ni, mi] = ind2sub([N,M], nmi);
    ppscore_single_row(1, nmi) = norm(baby_data(1,:).' - squeeze(DS_plane_feat(ni,mi,:)));
end

% No need to fill the first row of direction matrix as it already has
% zeros. 
%
% The first row of direction matrix is filled with 0is indicating this is the 
% starting point and that there are no previous points. 
% direction [1, :] = 0

% [Rules to fill the other rows of the table]
for pi = 2:P    
    
    ppscore_prev_row = ppscore_single_row;
    
    % no need for this step as the values will be overwritten anyway inside
    % the loop;
    % ppscore_single_row = 0;

    for nmi = 1:(N*M)

        [ni, mi] = ind2sub([N,M], nmi);
        
        % This stores the LN norm of this datapoint and DSP point. 
        % This will be added to the best optimal partial trjectory’s score 
        % where this point will be added as a current new member. 
        dist = norm(baby_data(pi, :).' - squeeze(DS_plane_feat(ni, mi, :))); 
               
 
        [bpt_score, previous_point_index] = best_parital_traj_score(nmi, radius, ppscore_prev_row, N, M);

        

        
        % [Done]: For tracing we only use the direction matrix, therefore
        % we don't need to compute scores for the whole ppscore matrix. 
        % We only need to store that of the previous one. This will save memory.
        
        % partial score of the partial trajeotry if this point is the last
        % one so far in the trajecotry
        ppscore_single_row(1, nmi) = dist + bpt_score;

                
        % Adding the index of the previous point which will be used to 
        % trace the optimal trajectory.

        direction(pi,nmi) = previous_point_index;
        
    end
    
end

projection_points_in_DSP = trace_the_optimum_traj_withnan(direction, ppscore_single_row, baby_data_feat_withnan);

disp('OP found!!');
end



function projection_points_in_DSP = trace_the_optimum_traj_withnan(direction, ppscore_single_row, baby_data_withnan)

P = size(direction,1);

P2 = size(baby_data_withnan,1);


projection_points_in_DSP = ones(P2,1)*nan;

% first get the last point of the best trajectory
[minval, minvalidx] = min(ppscore_single_row);


if isnan(minval)
    return;
end

% projection_points_in_DSP(end) = minvalidx;

% now tracing the tractectory
pi = P;
for pi2 = (P2):-1:1

    if ~all(isfinite(baby_data_withnan(pi2, :)))
        continue;
    end
    
    if pi == P
        projection_points_in_DSP(pi2) = minvalidx;
        last_projection = minvalidx;
        pi = pi -1;
        continue;
    end
    
    
    projection_points_in_DSP(pi2) = direction(pi+1, last_projection);
    
    last_projection = projection_points_in_DSP(pi2);
    
    pi = pi -1;
end
 
assert(pi == 0);

assert(pi2 == 1);

end


function [bpt_score, previous_point_index] = best_parital_traj_score(nmi, radius, ppscore_vec, N, M)

% [TODO] Write best_parital_traj_scorethis function
% This function takes the current  DSP index (nmi) and radius. 
% These two can be used to find the indices of neighbours within 
% the search radius. From these indices this function finds the 
% index with minimum ppscore in the previous row (or previous point
% in this partial trajectory). This function returns the ppscore of
% the minimum and its index.

% Getting the search range
[ni, mi] = ind2sub([N,M], nmi);


if ni-radius < 1
    n_st = 1;
else
    n_st = ni-radius;
end

if ni+radius > N
    n_ed = N;
else
    n_ed = ni+radius;
end


if mi-radius < 1
    m_st = 1;
else
    m_st = mi-radius;
end

if mi+radius > M
    m_ed = M;
else
    m_ed = mi+radius;
end
            

%convert ppscore_vec to DSP dimension mat
ppscore_vec_mat = reshape(ppscore_vec, [N,M]);

                                   
% getting the location of the last point of the best partial
% trajectory
myNewMat = ppscore_vec_mat(n_st:n_ed, m_st:m_ed);

[bpt_score, idx] = nanmin(myNewMat(:));            


[a,b]  = ind2sub(size(myNewMat), idx);

a = a + n_st - 1;
b = b + m_st - 1;

previous_point_index = sub2ind([N,M], a,b);

end

function filtered_baby_data = filtering_finite_rows(baby_data)
% keeping only those row with all finite values

[ri, ~] = find(~isfinite(baby_data));

baby_data(unique(ri), :) = [];

filtered_baby_data = baby_data; 

end
