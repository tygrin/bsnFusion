function [ Trans ] = add2TransBsn(intervals, Trans)
% add points to an existing transition list
% Input:
%   points      a matrix, each row is a point 
%   Trans       an existing transition list 
%               each row records a transition formated as following
%               [ start point, end point ]
%   Grid        an existing Grid created by the function createGrid()
%   isGrid      a boolean value to indicate whether discretizing 
%               the embedding space
%
% Author:   Zhifei Zhang
% E-mail:   zzhang61@vols.utk.edu
% Date:     July 20th, 2016

if nargin < 2
    error('Not enough input arguments!')
end

% approximate the points to the nearest grid cell
temp = intervals;

% append the new transitions to the end of the transition list
Trans = [ Trans ; temp ];



