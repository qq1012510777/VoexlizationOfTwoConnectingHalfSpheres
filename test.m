clc
clear all
close all

NumPntsInSingleVexel = 500;

currentPath = fileparts(mfilename('fullpath'));
addpath(genpath([currentPath, '/include']));

stl_water = stlread('water.stl'); % % % change this name if

points_water = stl_water.Points;
faces_water = stl_water.ConnectivityList;

clear stl_water

range_ = [-5e-4 5e-4];

RandomPnts = [];

RandomPnts(:, 1) = unifrnd(range_(1), range_(2), [NumPntsInSingleVexel, 1]);
RandomPnts(:, 2) = unifrnd(range_(1), range_(2), [NumPntsInSingleVexel, 1]);
RandomPnts(:, 3) = unifrnd(range_(1), range_(2), [NumPntsInSingleVexel, 1]);

%in = intriangulation(points_water,faces_water,RandomPnts);
in = inpolyhedron(faces_water,points_water,RandomPnts);

figure(1)
view(3)
patch('Vertices', points_water, 'Faces', faces_water, 'FaceVertexCData', ...
    zeros(size(faces_water, 1), 1), ...
    'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, ...
    'edgecolor', 'r'); hold on
as = find(in == 1);
ak = find(in == 0);
scatter3(RandomPnts(as, 1), RandomPnts(as, 2), RandomPnts(as, 3), ...
    'o', 'b', 'filled'); hold on

scatter3(RandomPnts(ak, 1), RandomPnts(ak, 2), RandomPnts(ak, 3), ...
    'o', 'k', 'filled'); hold on