%% My
%% add path
clc; clear; close all
addpath('data', 'DE', 'MGM', 'utilities')


%% import 

load_javaplex
import edu.stanford.math.plex4.*;

%% load data (two datasets are available in this demo)
datasetInd = 1; % <<--- please specify the dataset index, 1 or 2.
dataset = {
    'UCI_CharacterTrajectories' % #1
%    'MSR_Action3D'              % #2
    };
disp([ 'Processing the ' dataset{datasetInd}  ' dataset ' ])


%% load data (two datasets are available in this demo)


%switch datasetInd
%    case 1
        
    load(dataset{1})
    data = UCI_CharacterTrajectories.data;
    trueLabel = UCI_CharacterTrajectories.trueLabel;
    categories = UCI_CharacterTrajectories.categories;
    % load setting
    setting_UCI
    % split training and testing sets
    % the first half for training, and the rest for testing 
    trainInd = 1:1433;
    testInd = 1434:length(trueLabel);
% case 2
% otherwise
% error('Wrong dataset index!')
% end



%% create grid for each class
classLabels = unique(trueLabel);
n_class = length(classLabels);
n_dimSignal = size(data{1},1);
Trans = cell(n_class, 1);

TranStruct.data = UCI_CharacterTrajectories.data(trainInd);
TranStruct.label = UCI_CharacterTrajectories.trueLabel(trainInd);

%% training
startTime_train = tic;
for loop = 1:length(trainInd)
    if mod(loop, print_period) == 0
        fprintf('Trained %d / %d\n', loop, length(trainInd))
    end
    % extract data and label
    x = data{trainInd(loop)};
    % low-pass filter
    for i = 1:size(x, 1)
        x(i, :) = lowpassFilter(x(i,:), filter_param);
    end

    y = trueLabel(trainInd(loop));  
    % multi-dimensional delay embedding
    point_cloud = delayEmbedingND(x', DE_dim, DE_step, DE_slid);  
    % computing the persistence diagram of point cloud
    max_dimension = 2;
    max_filtration_value = 5;%??
    num_divisions = size(x, 1);
    
    stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);
    num_simplices = stream.getSize();
    
    % get persistence algorithm over Z/2Z
    persistence = api.Plex4.getModularSimplicialAlgorithm(max_dimension, 2);
    
    % compute the intervals
    intervals = persistence.computeIntervals(stream);
    
   % intervals_dim0=intervals.getIntervalsAtDimension(0);
   % intervals_dim1=intervals.getIntervalsAtDimension(1);
   TranStruct.intervals(loop) = intervals;
    % create the barcode plots
     %options.filename = 'ripsTorus';
     %options.max_filtration_value = max_filtration_value;
     %options.max_dimension = max_dimension - 1;
     %options.side_by_side = true;
     %plot_barcodes(intervals, options);
    % UCI_CharacterTrajectories.persistence = cell{}    
     %Trans{y} = add2TransBsn(intervals, Trans{y});
    
        
end
endTime_train = toc(startTime_train);

%% testing
startTime_test = tic;
dist = zeros(n_class, 1);
prediction = zeros(length(testInd), 1);

for loop = 1:length(testInd)
    if mod(loop, print_period) == 0
        fprintf('tested %d / %d\n', loop, length(testInd))
    end
    % extract data and label
    x = data{testInd(loop)};

    % low-pass filter
    for i = 1:size(x, 1)
        x(i, :) = lowpassFilter(x(i,:), filter_param);
    end

    y = trueLabel(testInd(loop));
    % multi-dimensional delay embedding
    point_cloud = delayEmbedingND(x', DE_dim, DE_step, DE_slid);
    % computing the persistence diagram of point cloud
    max_dimension = 2;
    max_filtration_value = 5;%??
    num_divisions = size(x, 1);
    
    stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);
    num_simplices = stream.getSize();
    
    % get persistence algorithm over Z/2Z
    persistence = api.Plex4.getModularSimplicialAlgorithm(max_dimension, 2);
    
    % compute the intervals
    intervals = persistence.computeIntervals(stream);
    
   % intervals_dim0=intervals.getIntervalsAtDimension(0);
   % intervals_dim1=intervals.getIntervalsAtDimension(1);
   
    % model matching
% 	for i = 1:n_class
%         dist(i) = HDist( point_cloud, Trans{i}, Grid{i}, ...
%             alpha, beta, isGrid );
% 	end 
%     [~, loc] = min(dist);
%     prediction(loop) = loc; 
%     

    for i = 1:n_class
        trainClass.intervals = TranStruct.intervals(TranStruct.label == i);
        dgmDist = PSdistComputation(intervals, trainClass.intervals);
        dist(i) = dgmDist;
    end
        [~, loc] = min(dist);
        prediction(loop) = loc; 
end

endTime_test = toc(startTime_test);



%% print training and testing time
fprintf('Training time: %.3fsec, %.3fsec per sample\n', ...
    endTime_train, endTime_train/length(trainInd))
fprintf('Testing time: %.3fsec, %.3fsec per sample\n', ...
    endTime_test, endTime_test/length(testInd))

%% plot confusion matrix
[CM, fig_handle] = confusionMatrix(trueLabel(testInd), prediction, categories);
fig = figure(fig_handle);
Accuracy = mean(trueLabel(testInd)==prediction);
fprintf('Accuracy = %.2f%%\n', Accuracy*100);
strAccu = sprintf('Confusion Matrix of the %s dataset, Accuracy %.2f%%', ...
    dataset{datasetInd}, Accuracy*100);
strAccu = strrep(strAccu, '_', '\_');
title(strAccu, 'fontsize', 16);
set(gcf, 'units','normalized','outerposition',[.1 .1 .8 .8])


