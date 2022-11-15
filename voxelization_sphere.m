clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));
addpath(genpath([currentPath, '/include']));

NumPntsInSingleVexel = 100;
criterionVolume = 0.5;
gridsize = 2.5e-4; % side length of a voxel
NUMthreads = 10; % number of threads

% load the surface mesh of the model
% load the surface mesh of the model
% load the surface mesh of the model
stl = stlread('model.stl'); % % % change this name

points = stl.Points;
faces = stl.ConnectivityList;
clear stl;

minX = min(points(:, 1));
maxX = max(points(:, 1));
minY = min(points(:, 2));
maxY = max(points(:, 2));
minZ = min(points(:, 3));
maxZ = max(points(:, 3));

if (minX ~= minY || minX ~= minZ || maxX ~= maxY || maxX ~= maxZ)
    disp([minX, maxX, minY, maxY, minZ, maxZ])
    error('The initial voxelization should be cubic')
end

% voexlization of a large cubic domain that encloses the model
% voexlization of a large cubic domain that encloses the model
% voexlization of a large cubic domain that encloses the model
% the voexlization consists of many cubes
% the voexlization consists of many cubes
% the voexlization consists of many cubes
[X, Y, Z] = meshgrid([minX:gridsize:maxX], [minY:gridsize:maxY], [minZ:gridsize:maxZ]);

NUMVertices = size(X(:), 1);

[structure___cube] = CreateStructureCubes(X, Y, Z);
Vpoints = [X(:), Y(:), Z(:)]; % points of all cubes' vertices
clear X Y Z

figure(1);
view(3);
title('Voxelization of a big, cubic domain enclosing the model', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(:, [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0); hold on;

% Vpoints_1_ = Vpoints(find(Vpoints(:, 3) > 2.07e-19), :);
% Vpoints_2_ = Vpoints(find(Vpoints(:, 3) < 2.09e-19), :);

% identify the upper and lower half spheres
% identify the upper and lower half spheres
% identify the upper and lower half spheres
ks = points(:, 3);

% because the connection point between the upper and lower half spheres is
% 0, 0, 0
ks(find(abs(ks) < 1e-7), 1) = 0;

ks(find(ks > 0), 1) = 1;
ks(find(ks < 0), 1) = 2;

Face_att(:, 1) = ks(faces(:, 1), 1);
Face_att(:, 2) = ks(faces(:, 2), 1);
Face_att(:, 3) = ks(faces(:, 3), 1);

Face_att = Face_att(:, 1);

% the following is the faces (triangles) of
% upper and lower half spheres
Face1_1_ = find(Face_att == 1);
Face1_2_ = find(Face_att == 2);

% now the convex hull of upper half sphere
fs_1 = faces(Face1_1_, :);
convex_hull_1 = points (fs_1(:), :);

% now the convex hull of lower half sphere
fs_2 = faces(Face1_2_, :);
convex_hull_2 = points (fs_2(:), :);

% identify if a voxel point belong to the half spheres
% identify if a voxel point belong to the half spheres
% identify if a voxel point belong to the half spheres
b_1_ = in_convex_polyhedron(convex_hull_1, Vpoints([1:NUMVertices], :), zeros(NUMVertices, 1), NUMthreads);
s_1_ = find(b_1_ == 1);

b_2_ = in_convex_polyhedron(convex_hull_2, Vpoints([1:NUMVertices], :), b_1_, NUMthreads);
s_2_ = find(b_2_ == 1);

figure(2)
title('Points which are inside the upper, half sphere', 'interpreter', 'latex')
patch('Vertices', points, 'Faces', faces(Face1_1_, :), 'FaceVertexCData', points(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'facealpha', 0); view(3); colorbar; hold on;
scatter3(Vpoints(s_1_, 1), Vpoints(s_1_, 2), Vpoints(s_1_, 3), '.'); hold on
figure(3)
title('Points which are inside the lower, half sphere', 'interpreter', 'latex')
patch('Vertices', points, 'Faces', faces(Face1_2_, :), 'FaceVertexCData', points(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'facealpha', 0); view(3); colorbar; hold on;
scatter3(Vpoints(s_2_, 1), Vpoints(s_2_, 2), Vpoints(s_2_, 3), '.'); hold on

% now let check if one voxel is inside the spheres
% now let check if one voxel is inside the spheres
% now let check if one voxel is inside the spheres

cube_1 = [];
cube_2 = [];

cube_1_i = [];
cube_2_i = [];

VerticesID = [1:8];

for i = 1:size(structure___cube, 1)
    display(['check cube NO ', num2str(i), '/', num2str(size(structure___cube, 1))]);

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

    k_1_ = in_convex_polyhedron(convex_hull_1, RandomPnts, zeros(NumPntsInSingleVexel, 1), NUMthreads);
    e_1_ = find(k_1_ == 1);

    VolumeFraction1 = size(e_1_, 1) / size(k_1_, 1);
    % close a7
    if (VolumeFraction1 >= criterionVolume)
        cube_1_i = [cube_1_i; i];
    end

    % monte carlo method to identify the intersection volume between cube and
    % upper half spehere, if just some of vertices, or even no vertices,
    % are inside the lower half sphere
    k_2_ = in_convex_polyhedron(convex_hull_2, RandomPnts, zeros(NumPntsInSingleVexel, 1), NUMthreads);
    e_2_ = find(k_2_ == 1);

    VolumeFraction2 = size(e_2_, 1) / size(k_2_, 1);

    if (VolumeFraction2 >= criterionVolume)
        cube_2_i = [cube_2_i; i];
    end

    clear minX_t maxX_t minY_t maxY_t minZ_t maxZ_t RandomPnts k_1_ e_1_ k_2_ e_2_
end

% load('resolution_5e-5_ff.mat')

figure(4);
view(3);
xlabel('x')
ylabel('y')
zlabel('z'); hold on
title('cubes which are completely inside the upper and lower, half sphere', 'interpreter', 'latex')
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

figure(5);
view(3);
xlabel('x')
ylabel('y')
zlabel('z'); hold on
title('cubes which are partialy inside the upper, half sphere', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_1_i, [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'r'); hold on;

figure(6);
view(3);
xlabel('x')
ylabel('y')
zlabel('z'); hold on
title('cubes which are partialy inside the lower, half sphere', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube(cube_2_i, [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0.0, 'edgecolor', 'b'); hold on;

figure(7);
view(3);
xlabel('x')
ylabel('y')
zlabel('z'); hold on
title('The voxelization of the input model', 'interpreter', 'latex')
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [1:4]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [5:8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [1, 2, 6, 5]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [2, 3, 7, 6]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [3, 4, 8, 7]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;
patch('Vertices', Vpoints, 'Faces', structure___cube([cube_1; cube_2; cube_1_i; cube_2_i], [4, 1, 5, 8]), 'FaceVertexCData', Vpoints(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1.0); hold on;

% write nas file
% write nas file
% write nas file
precision_string = '%6.4f';
fid = fopen([currentPath, '/HexahedronMesh.nas'], 'w');
fprintf(fid, "$ Generated by TC Yin\n");
fprintf(fid, "$ ");
fprintf(fid, [date, '\n']);
fprintf(fid, "BEGIN BULK\n");
fprintf(fid, "$ Grid data section\n");

for i = 1:NUMVertices %size(Vpoints, 1)

    fprintf(fid, '%-8s%-8s%-8s%8s%8s%8s\n', 'GRID', num2str(i, 6), '', ...
        ... %num2str(Vpoints(i, 1), precision_string), num2str(Vpoints(i, 2), precision_string), num2str(Vpoints(i, 3), precision_string));
        Num2Str_Set8Width(Vpoints(i, 1)), Num2Str_Set8Width(Vpoints(i, 2)), Num2Str_Set8Width(Vpoints(i, 3)));
end

fprintf(fid, "$ Element data section\n");

gy = 1;

for i = 1:size(structure___cube, 1)

    a = find([cube_1; cube_2; cube_1_i; cube_2_i] == i);

    if (isempty(a) == 0)
        fprintf(fid, '%-8s%-8s%-8s%8s%8s%8s%8s%8s%8s+CONT\n', 'CHEXA', num2str(gy, 5), '1', ...
            num2str(structure___cube(i, 1), '%6d'), ...
            num2str(structure___cube(i, 4), '%6d'), ...
            num2str(structure___cube(i, 3), '%6d'), ...
            num2str(structure___cube(i, 2), '%6d'), ...
            num2str(structure___cube(i, 5), '%6d'), ...
            num2str(structure___cube(i, 8), '%6d'));
        fprintf(fid, '%-8s%-8s%-8s\n', '+CONT', ...
            num2str(structure___cube(i, 7), '%6d'), ...
            num2str(structure___cube(i, 6), '%6d'));
        gy = gy + 1;
    end

end

fprintf(fid, "ENDDATA\n");
fclose(fid);
