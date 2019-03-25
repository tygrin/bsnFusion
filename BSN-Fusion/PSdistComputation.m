function [ Dist0] = PSdistComputation(intervals, intervals_Train)
% compute the ps distance between the given trajectory 
% (points) and a learned model (Trans).
% Input:
%   interval      test sample intervel 
%   intervals_Train      training sample intervals

 m = length(intervals_Train); 
 
 dist_dim0 = zeros(m, 1);
 %dist_dim1 = zeros(m, 1);
 
 for i=1:m
     intervalsA = intervals;
     intervalsB = intervals_Train(i);
     intervalsA_dim0=intervalsA.getIntervalsAtDimension(0);
     intervalsB_dim0=intervalsB.getIntervalsAtDimension(0);
     bottleneck_distance_dim0 = edu.stanford.math.plex4.bottleneck.BottleneckDistance.computeBottleneckDistance(intervalsA_dim0,intervalsB_dim0);
     
     dist_dim0(i) = bottleneck_distance_dim0;
     
     %intervalsA_dim1=intervalsA.getIntervalsAtDimension(1);
     %intervalsB_dim1=intervalsB.getIntervalsAtDimension(1);
     %bottleneck_distance_dim1 = edu.stanford.math.plex4.bottleneck.BottleneckDistance.computeBottleneckDistance(intervalsA_dim1,intervalsB_dim1);
     
     %dist_dim1(i) = bottleneck_distance_dim1;
 end
 
    Dist0 = mean(dist_dim0);
    %Dist1 = min(dist_dim1);
 
% Dim0 computation
 
% [p, ~] = size(Trans);
% 

% if nargin < 6
%     isGrid = true;
% end
% if nargin < 5
%     beta = 1.0;
% end
% if nargin < 4
%     alpha = 1.0;
% end
% if nargin < 3
%     error('Not enough input arguments!')
% end
% if isempty(Trans) || isempty(Grid)
%     error('Transition list or Grid is empty!')
% end
% 

% % approximate the points to the nearest grid cell
% if isGrid
%     gridCenter = Grid.center;
%     gridSize = Grid.size;
%     points = round((points-repmat(gridCenter,m,1)) ./ ...
%         repmat(gridSize,m,1)) .* repmat(gridSize,m,1) + ...
%         repmat(gridCenter,m,1);
% end
% 
% % direction, location and length of transitions in the embedding space
% vec_Trans = Trans(:, n+1:2*n) - Trans(:, 1:n);
% loc_Trans = (Trans(:, n+1:2*n) + Trans(:, 1:n)) / 2;
% len_Trans = sqrt(sum(vec_Trans.^2, 2));
% 
% % direction, location and length of the given trajectory
% vec_points = points(2:end, :) - points(1:end-1, :);
% loc_points = (points(2:end, :) + points(1:end-1, :)) / 2;
% len_points = sqrt(sum(vec_points.^2, 2));
% 
% % normalized angle between learned transitions and given trajectory
% norm_angle = exp( real(acos(vec_points*vec_Trans' ./ ...
%     (len_points*len_Trans') )));
% norm_angle(len_points==0, len_Trans==0) = 0;
% 
% % normalized lenght difference
% norm_length = exp( (repmat(len_points, 1, p) - ...
%     repmat(len_Trans', m-1, 1)).^2 ./ (repmat(len_points, 1, p).^2) );
% norm_length(isnan(norm_length)) = 0;
% 
% % modified Hausdorff distance
% norm_distance = zeros(m-1, p);
% for i=1:m-1
%     if len_points(i)>0
%         norm_distance(i, :) = sqrt(sum((repmat(loc_points(i,:), p, 1)...
%             - loc_Trans).^2, 2))' / len_points(i);
%     else
%         norm_distance(i, :) = sqrt(sum((repmat(loc_points(i,:), p, 1)...
%             - loc_Trans).^2, 2))';
%     end
% end
% norm_dist = norm_distance + alpha*norm_length + beta*norm_angle;
% dist = min(norm_dist, [], 2);
% dist = mean( dist(len_points > 0) );


