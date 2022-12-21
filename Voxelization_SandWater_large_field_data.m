clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));
addpath(genpath([currentPath, '/include']));

NumPntsInSingleVexel = 200; % number of points in a 3D voxel
criterionVolume = 0.5; % the critical fraction if a voxel is belong to the object
gridsize = 1e-4; % side length of a voxel
% NUMthreads = 10; % number of threads
offset_vertical_fraction_z = 0.3; % offset value (in form of percentage), for example 0.5 = 50 % of gridsize
offset_vertical_fraction_y = 0.4; % offset value (in form of percentage), for example 0.5 = 50 % of gridsize
offset_vertical_fraction_x = 0.5; % offset value (in form of percentage), for example 0.5 = 50 % of gridsize

if ((offset_vertical_fraction_z < 0 || offset_vertical_fraction_z >= 1) || ...
        (offset_vertical_fraction_y < 0 || offset_vertical_fraction_y >= 1) || ...
        (offset_vertical_fraction_x < 0 || offset_vertical_fraction_x >= 1))
    error('Wrong offset value!');
end

stl_sand = stlread('sand.stl'); % % % change this name if
stl_water = stlread('water.stl'); % % % change this name if

points_sand = stl_sand.Points;
faces_sand = stl_sand.ConnectivityList;

points_water = stl_water.Points;
faces_water = stl_water.ConnectivityList;

clear stl_sand  stl_water

figure(1)
view(3)
patch('Vertices', points_sand, 'Faces', faces_sand, 'FaceVertexCData', zeros(size(faces_sand, 1), 1), ...
    'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, ...
    'edgecolor', 'r'); hold on
pbaspect([1, 1, 1]); hold on
xlabel('x')
ylabel('y')
zlabel('z'); hold on

figure(2)
view(3)
patch('Vertices', points_water, 'Faces', faces_water, 'FaceVertexCData', zeros(size(faces_water, 1), 1), ...
    'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, ...
    'edgecolor', 'b'); hold on
pbaspect([1, 1, 1]); hold on
xlabel('x')
ylabel('y')
zlabel('z'); hold on

% let's address the sand first
% let's address the sand first
% let's address the sand first
minX = min(points_sand(:, 1));
maxX = max(points_sand(:, 1));
minY = min(points_sand(:, 2));
maxY = max(points_sand(:, 2));
minZ = min(points_sand(:, 3));
maxZ = max(points_sand(:, 3));

if (minX ~= minY || minX ~= minZ || maxX ~= maxY || maxX ~= maxZ)
    disp([minX, maxX, minY, maxY, minZ, maxZ])
    error('The initial voxelization should be cubic')
end

% mesh

X = []; Y = []; Z = [];
Ax = []; Ay = []; Az = [];

if (offset_vertical_fraction_x ~= 0)
    Ax = [minX, (minX + (1 - offset_vertical_fraction_x) * gridsize):gridsize:(maxX - offset_vertical_fraction_x * gridsize), ...
            maxX];
else
    Ax = [minX:gridsize:maxX];
end

if (offset_vertical_fraction_y ~= 0)
    Ay = [minY, (minY + (1 - offset_vertical_fraction_y) * gridsize):gridsize:(maxY - offset_vertical_fraction_y * gridsize), ...
            maxX];
else
    Ay = [minY:gridsize:maxY];
end

if (offset_vertical_fraction_z ~= 0)
    Az = [minZ, (minZ + (1 - offset_vertical_fraction_z) * gridsize):gridsize:(maxZ - offset_vertical_fraction_z * gridsize), ...
            maxZ];
else
    Az = [minZ:gridsize:maxZ];
end

[X, Y, Z] = meshgrid(Ax, Ay, Az);

NUMVertices = size(X(:), 1);

[structure___cube] = CreateStructureCubes(X, Y, Z);
Vpoints = [X(:), Y(:), Z(:)]; % points of all cubes' vertices
clear X Y Z

figure(3);
view(3);
title('Voxelization of a big, cubic domain enclosing the model (sand)', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
pbaspect([1, 1, 1]); hold on
xlabel('x')
ylabel('y')
zlabel('z'); hold on
% identify the upper and lower half spheres
% identify the upper and lower half spheres
% identify the upper and lower half spheres
ks = points_sand(:, 3);

% because the connection point between the upper and lower half spheres is
% 0, 0, 0
ks(find(abs(ks) < 1e-7), 1) = 0;

ks(find(ks > 0), 1) = 1;
ks(find(ks < 0), 1) = 2;

Face_att(:, 1) = ks(faces_sand(:, 1), 1);
Face_att(:, 2) = ks(faces_sand(:, 2), 1);
Face_att(:, 3) = ks(faces_sand(:, 3), 1);

Face_att = Face_att(:, 1);

% the following is the faces (triangles) of
% upper and lower half spheres
Face1_1_ = find(Face_att == 1);
Face1_2_ = find(Face_att == 2);

% now the convex hull of upper half sphere
fs_1 = faces_sand(Face1_1_, :);
% convex_hull_1 = points_sand (fs_1(:), :);

% now the convex hull of lower half sphere
fs_2 = faces_sand(Face1_2_, :);
% convex_hull_2 = points_sand (fs_2(:), :);

% identify if a voxel point belong to the half spheres
% identify if a voxel point belong to the half spheres
% identify if a voxel point belong to the half spheres
% b_1_ = in_convex_polyhedron(convex_hull_1, Vpoints([1:NUMVertices], :), zeros(NUMVertices, 1), NUMthreads);

b_1_ = intriangulation(points_sand, faces_sand(Face1_1_, :), Vpoints([1:NUMVertices], :));
s_1_ = find(b_1_ == 1);

%b_2_ = in_convex_polyhedron(convex_hull_2, Vpoints([1:NUMVertices], :), b_1_, NUMthreads);
b_2_ = intriangulation(points_sand, faces_sand(Face1_2_, :), Vpoints([1:NUMVertices], :));
s_2_ = find(b_2_ == 1);

figure(4)
title('Points which are inside the upper, half sphere (sand)', 'interpreter', 'latex')
patch('Vertices', points_sand, 'Faces', faces_sand(Face1_1_, :), 'FaceVertexCData', points_sand(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'facealpha', 0); view(3); colorbar; hold on;
scatter3(Vpoints(s_1_, 1), Vpoints(s_1_, 2), Vpoints(s_1_, 3), '.'); hold on
pbaspect([1, 1, 1]); hold on
xlabel('x')
ylabel('y')
zlabel('z'); hold on
figure(5)
title('Points which are inside the lower, half sphere (sand)', 'interpreter', 'latex')
patch('Vertices', points_sand, 'Faces', faces_sand(Face1_2_, :), 'FaceVertexCData', points_sand(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'facealpha', 0); view(3); colorbar; hold on;
scatter3(Vpoints(s_2_, 1), Vpoints(s_2_, 2), Vpoints(s_2_, 3), '.'); hold on
pbaspect([1, 1, 1]); hold on
xlabel('x')
ylabel('y')
zlabel('z'); hold on

% now let check if one voxel is inside the spheres
% now let check if one voxel is inside the spheres
% now let check if one voxel is inside the spheres

cube_1 = [];
cube_2 = [];

cube_1_i = [];
cube_2_i = [];

VerticesID = [1:8];

for i = 1:size(structure___cube, 1)
    display(['check (sand) cube NO ', num2str(i), '/', num2str(size(structure___cube, 1))]);

    % if all the vertices of a cube is completely inside the upper half sphere
    a = ismember(structure___cube(i, VerticesID), s_1_);

    if (sum(a) == size(a, 2))
        cube_1 = [cube_1; i];
        continue
    end

    % if all the vertices of a cube is completely inside the lower half sphere
    b = ismember(structure___cube(i, VerticesID), s_2_);

    if (sum(b) == size(b, 2))
        cube_2 = [cube_2; i];
        continue
    end

    % monte carlo method to identify the intersection volume between cube and
    % upper half spehere, if just some of vertices, or even no vertices,
    % are inside the upper half sphere
    RandomPnts = zeros(NumPntsInSingleVexel, 3);
    minX_t = Vpoints(structure___cube(i, 1), 1);
    maxX_t = Vpoints(structure___cube(i, 4), 1);
    RandomPnts(:, 1) = unifrnd(minX_t, maxX_t, [NumPntsInSingleVexel, 1]);
    %
    minY_t = Vpoints(structure___cube(i, 1), 2);
    maxY_t = Vpoints(structure___cube(i, 2), 2);
    RandomPnts(:, 2) = unifrnd(minY_t, maxY_t, [NumPntsInSingleVexel, 1]);
    %
    minZ_t = Vpoints(structure___cube(i, 1), 3);
    maxZ_t = Vpoints(structure___cube(i, 5), 3);
    RandomPnts(:, 3) = unifrnd(minZ_t, maxZ_t, [NumPntsInSingleVexel, 1]);

    if (minX_t >= maxX_t)
        error(['Incorrect voxel range:\n', num2str([minX_t miaxX_t]), '\n', num2str([minY_t miaxY_t]), '\n', num2str([minZ_t miaxZ_t]), '\n'])
    end

    %     figure(a7);
    %     view(3);
    %     patch('Vertices', Vpoints, 'Faces', structure___cube(i, [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
    %     patch('Vertices', Vpoints, 'Faces', structure___cube(i, [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
    %     patch('Vertices', Vpoints, 'Faces', structure___cube(i, [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
    %     patch('Vertices', Vpoints, 'Faces', structure___cube(i, [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
    %     patch('Vertices', Vpoints, 'Faces', structure___cube(i, [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
    %     patch('Vertices', Vpoints, 'Faces', structure___cube(i, [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
    %     scatter3(RandomPnts(:, 1), RandomPnts(:, 2), RandomPnts(:, 3), 'k', 'filled'); hold on
    %     pbaspect([1, 1, 1]); hold on
    %     xlabel('x')
    %     ylabel('y')
    %     zlabel('z'); hold on

    %k_1_ = in_convex_polyhedron(convex_hull_1, RandomPnts, zeros(NumPntsInSingleVexel, 1), NUMthreads);
    k_1_ = intriangulation(points_sand, faces_sand(Face1_1_, :), RandomPnts);
    e_1_ = find(k_1_ == 1);

    VolumeFraction1 = size(e_1_, 1) / size(k_1_, 1);
    % close a7
    if (VolumeFraction1 >= criterionVolume)
        cube_1_i = [cube_1_i; i];
    end

    % monte carlo method to identify the intersection volume between cube and
    % upper half spehere, if just some of vertices, or even no vertices,
    % are inside the lower half sphere
    % k_2_ = in_convex_polyhedron(convex_hull_2, RandomPnts, zeros(NumPntsInSingleVexel, 1), NUMthreads);
    k_2_ = intriangulation(points_sand, faces_sand(Face1_2_, :), RandomPnts);
    e_2_ = find(k_2_ == 1);

    VolumeFraction2 = size(e_2_, 1) / size(k_2_, 1);

    if (VolumeFraction2 >= criterionVolume)
        cube_2_i = [cube_2_i; i];
    end

    clear minX_t maxX_t minY_t maxY_t minZ_t maxZ_t RandomPnts k_1_ e_1_ k_2_ e_2_
end

% load('resolution_5e-5_ff.mat')

figure(6);
view(3);
xlabel('x')
ylabel('y')
zlabel('z'); hold on
pbaspect([1, 1, 1]); hold on
title('cubes which are completely inside the upper and lower, half sphere (sand)', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1, [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1, [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1, [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1, [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1, [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1, [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;

patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2, [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2, [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2, [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2, [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2, [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2, [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;

figure(7);
view(3);
xlabel('x')
ylabel('y')
zlabel('z'); hold on
pbaspect([1, 1, 1]); hold on
title('cubes which are partialy inside the upper, half sphere (sand)', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;

figure(8);
view(3);
xlabel('x')
ylabel('y')
zlabel('z'); hold on
pbaspect([1, 1, 1]); hold on
title('cubes which are partialy inside the lower, half sphere (sand)', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;

figure(9);
view(3);
xlabel('x')
ylabel('y')
zlabel('z'); hold on
pbaspect([1, 1, 1]); hold on
title('The voxelization of the input model (sand)', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;

%----------------
% let's address water
%----------------

cube_3 = [1:1:size(structure___cube, 1)]';
cube_3([cube_1; cube_2; cube_1_i; cube_2_i], :) = [];
cube_3_i = []; % for water cubes

minX_water = min(points_water(:, 1));
maxX_water = max(points_water(:, 1));

minY_water = min(points_water(:, 2));
maxY_water = max(points_water(:, 2));

minZ_water = min(points_water(:, 3));
maxZ_water = max(points_water(:, 3));

for i = 1:size(cube_3, 1)
    display(['check (water) cube NO ', num2str(i), '/', num2str(size(cube_3, 1))]);

    cubeNO_t = cube_3(i);

    % monte carlo method to identify the intersection volume between cube and
    % water, if just some of vertices, or even no vertices,
    % are inside the water stl
    RandomPnts = zeros(NumPntsInSingleVexel, 3);
    minX_t = Vpoints(structure___cube(cubeNO_t, 1), 1);
    maxX_t = Vpoints(structure___cube(cubeNO_t, 4), 1);
    RandomPnts(:, 1) = unifrnd(minX_t, maxX_t, [NumPntsInSingleVexel, 1]);
    %
    minY_t = Vpoints(structure___cube(cubeNO_t, 1), 2);
    maxY_t = Vpoints(structure___cube(cubeNO_t, 2), 2);
    RandomPnts(:, 2) = unifrnd(minY_t, maxY_t, [NumPntsInSingleVexel, 1]);
    %
    minZ_t = Vpoints(structure___cube(cubeNO_t, 1), 3);
    maxZ_t = Vpoints(structure___cube(cubeNO_t, 5), 3);
    RandomPnts(:, 3) = unifrnd(minZ_t, maxZ_t, [NumPntsInSingleVexel, 1]);

    if ((minZ_t < minZ_water && maxZ_t <= minZ_water) || ...
            (minZ_t >= maxZ_water && maxZ_t > maxZ_water))
        continue
    end

    if ((minY_t < minY_water && maxY_t <= minY_water) || ...
            (minY_t >= maxY_water && maxY_t > maxY_water))
        continue
    end

    if ((minX_t < minX_water && maxX_t <= minX_water) || ...
            (minX_t >= maxX_water && maxX_t > maxX_water))
        continue
    end

    if (minX_t >= maxX_t)
        error(['Incorrect voxel range:\n', num2str([minX_t miaxX_t]), '\n', num2str([minY_t miaxY_t]), '\n', num2str([minZ_t miaxZ_t]), '\n'])
    end

    in = inpolyhedron(faces_water, points_water, RandomPnts);

    e_1_ = find(in == 1);

    VolumeFraction1 = size(e_1_, 1) / NumPntsInSingleVexel;

    if (VolumeFraction1 >= criterionVolume)
        cube_3_i = [cube_3_i; cubeNO_t];

        %         figure(13);
        %         view(3);
        %         patch('Vertices', Vpoints, 'Faces', structure___cube(cubeNO_t, [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
        %         patch('Vertices', Vpoints, 'Faces', structure___cube(cubeNO_t, [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
        %         patch('Vertices', Vpoints, 'Faces', structure___cube(cubeNO_t, [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
        %         patch('Vertices', Vpoints, 'Faces', structure___cube(cubeNO_t, [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
        %         patch('Vertices', Vpoints, 'Faces', structure___cube(cubeNO_t, [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
        %         patch('Vertices', Vpoints, 'Faces', structure___cube(cubeNO_t, [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
        %         scatter3(RandomPnts(e_1_, 1), RandomPnts(e_1_, 2), RandomPnts(e_1_, 3), 'k', 'filled'); hold on
        %         e_2_ = find(in == 0);
        %         scatter3(RandomPnts(e_2_, 1), RandomPnts(e_2_, 2), RandomPnts(e_2_, 3), 'r', 'filled'); hold on
        %         pbaspect([1 1 1])
        %         hold on
        %         patch('Vertices', points_water, 'Faces', faces_water, 'FaceVertexCData', zeros(size(faces_water, 1), 1), ...
        %             'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, ...
        %             'edgecolor', 'b'); hold on
        %         close 13
    end

end

figure(10);
view(3);
xlabel('x')
ylabel('y')
zlabel('z'); hold on
pbaspect([1, 1, 1]); hold on
title('The voxelization of water (water)', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
hold on
patch('Vertices', points_water, 'Faces', faces_water, 'FaceVertexCData', ...
    zeros(size(faces_water, 1), 1), ...
    'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, ...
    'edgecolor', 'b'); hold on

figure(11);
view(3);
xlabel('x')
ylabel('y')
zlabel('z'); hold on
pbaspect([1, 1, 1]); hold on
title('The voxelization of sand and water', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;

hold on
facecolor_alpha = 0.2;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'r', 'EdgeAlpha', 1, 'facealpha', facecolor_alpha); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'r', 'EdgeAlpha', 1, 'facealpha', facecolor_alpha); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'r', 'EdgeAlpha', 1, 'facealpha', facecolor_alpha); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'r', 'EdgeAlpha', 1, 'facealpha', facecolor_alpha); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'r', 'EdgeAlpha', 1, 'facealpha', facecolor_alpha); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_3_i], [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'r', 'EdgeAlpha', 1, 'facealpha', facecolor_alpha); hold on;

% write nas file
% write nas file
% write nas file
% load([currentPath, '/Voxelization_SandWater_16_Nov_2022_19_14_26_tt.mat'])
% currentPath = fileparts(mfilename('fullpath'));
% precision_string = '%6.4f';
fid = fopen([currentPath, '/HexahedronMeshSandAndWater.nas'], 'w');
fprintf(fid, "$ Generated by TC Yin\n");
fprintf(fid, "$ ");
fprintf(fid, [date, '\n']);
fprintf(fid, "BEGIN BULK\n");
fprintf(fid, "$ Grid data section\n");

Data_num_significant = 8;

for i = 1:NUMVertices %size(Vpoints, 1)
    disp(['write points: ', num2str(i), '/', num2str(NUMVertices)]);

    fprintf(fid, '%-8s%-16s%-16s%16s%16s%8s\n', 'GRID*', num2str(i, Data_num_significant), '', ...
        Num2Str_Set_Width(Vpoints(i, 1), Data_num_significant), Num2Str_Set_Width(Vpoints(i, 2), Data_num_significant), '*GRID');

    fprintf(fid, '%-8s%-16s\n', '*GRID', Num2Str_Set_Width(Vpoints(i, 3), Data_num_significant));

end

fprintf(fid, "$ Element data section\n");

gy = 1;

sandcubeNO = [cube_1; cube_2; cube_1_i; cube_2_i];

for j = 1:size(sandcubeNO, 1)
    i = sandcubeNO(j);
    disp(['write sand cubes: ', num2str(j), '/', num2str(size(sandcubeNO, 1))]);

    fprintf(fid, '%-8s%-16s%-16s%16s%16s%16s+CONT\n', 'CHEXA*', ...
        num2str(gy, Data_num_significant), ...
        '1', ...
        num2str(structure___cube(i, 1), '%16d'), ...
        num2str(structure___cube(i, 4), '%16d'), ...
        num2str(structure___cube(i, 3), '%16d'));
    fprintf(fid, '%-8s%16s%16s%16s%16s%16s\n', '+CONT', ...
        num2str(structure___cube(i, 2), '%16d'), ...
        num2str(structure___cube(i, 5), '%16d'), ...
        num2str(structure___cube(i, 8), '%16d'), ...
        num2str(structure___cube(i, 7), '%16d'), ...
        num2str(structure___cube(i, 6), '%16d'));
    gy = gy + 1;

end

for j = 1:size(cube_3_i, 1)
    disp(['write water cubes: ', num2str(j), '/', num2str(size(cube_3_i, 1))]);
    i = cube_3_i(j);

    fprintf(fid, '%-8s%-16s%-16s%16s%16s%16s+CONT\n', 'CHEXA*', ...
        num2str(gy, Data_num_significant), ...
        '2', ...
        num2str(structure___cube(i, 1), '%16d'), ...
        num2str(structure___cube(i, 4), '%16d'), ...
        num2str(structure___cube(i, 3), '%16d'));
    fprintf(fid, '%-8s%16s%16s%16s%16s%16s\n', '+CONT', ...
        num2str(structure___cube(i, 2), '%16d'), ...
        num2str(structure___cube(i, 5), '%16d'), ...
        num2str(structure___cube(i, 8), '%16d'), ...
        num2str(structure___cube(i, 7), '%16d'), ...
        num2str(structure___cube(i, 6), '%16d'));
    gy = gy + 1;
end

fprintf(fid, "ENDDATA\n");
fclose(fid);

sg = datestr(datetime);
asf = find(sg == ' ');
sg(asf) = '_';
asf = find(sg == '-');
sg(asf) = '_';
asf = find(sg == ':');
sg(asf) = '_';

save([currentPath, '/Voxelization_SandWater_', sg, '_tt.mat'])
