function [structure_] = CreateStructureCubes(X_1, Y_1, Z_1)
    % cubes
    NumCubes = size(X_1, 1) - 1; NumCubes = NumCubes.^3;

    sizeX = size(X_1, 1) - 1;

    structure_lowest = [];

    % the lowest mesh
    for i = 1:NumCubes / sizeX

        j = floor((i - 1) * 1.0 / sizeX) + i;

        structure_lowest = [structure_lowest; j, j + 1, (j + 2) + sizeX, (j + 1) + sizeX];

    end

    structure_ = zeros(NumCubes, 8);

    for i = 1:sizeX

        structure_([(i - 1) * (NumCubes / sizeX) + 1:i * (NumCubes / sizeX)], [1:4]) = ...
            structure_lowest + (i - 1) * size(X_1, 1) * size(X_1, 1);
        structure_([(i - 1) * (NumCubes / sizeX) + 1:i * (NumCubes / sizeX)], [5:8]) = ...
            structure_lowest + (i) * size(X_1, 1) * size(X_1, 1);

    end

%     % add mid points to edges
%     structure_tmp = zeros(size(structure_, 1), 20);
%     structure_tmp(:, [1, 3, 5, 7, 9, 11, 13, 15]) = structure_;
%     clear structure_
%     structure_ = structure_tmp;
%     clear structure_tmp
% 
%     points_new = [NaN, 0, 0];
% 
%     ksd = [1, 3, 5, 7];
%     ksf = [9, 11, 13, 15];
%     X = X_1(:);
%     Y = Y_1(:);
%     Z = Z_1(:);
%     
%     Vpoints = [X, Y, Z];
%     sizeOfPnts_now = size(Vpoints, 1);
%     for i = 1:size(structure_, 1)
%         disp(['add midpoints ', num2str(i), '/', num2str(size(structure_, 1)), ' ...'])
%         % --------------------------- top
%         for j = 1:4
%             node1 = structure_(i, ksd(j));
%             node2 = structure_(i, ksd(mod(j, 4) + 1));
% 
%             midpoint = 0.5 * [X(node1, 1) + X(node2, 1), ...
%                             Y(node1, 1) + Y(node2, 1), ...
%                             Z(node1, 1) + Z(node2, 1)];
%             isPresent = find(ismember(points_new, midpoint, 'rows'));
% 
%             if (isempty(isPresent))
%                 if(isnan(points_new(1, 1)) == 0)
%                     points_new = [points_new; midpoint];
%                 else
%                     points_new = midpoint;
%                 end
%                 structure_(i, ksd(j) + 1) = size(points_new, 1) + sizeOfPnts_now;
%             else
% 
%                 if (size(isPresent, 1) == 1)
%                     structure_(i, ksd(j) + 1) = isPresent(1, 1) + sizeOfPnts_now;
%                 else
%                     error('Duplicate points exist in array of points_new 1')
%                 end
% 
%             end
% 
%             clear midpoint
%         end
% 
%         % --------------------------- bottom
%         for j = 1:4
%             node1 = structure_(i, ksf(j));
%             node2 = structure_(i, ksf(mod(j, 4) + 1));
% 
%             midpoint = 0.5 * [X(node1, 1) + X(node2, 1), ...
%                             Y(node1, 1) + Y(node2, 1), ...
%                             Z(node1, 1) + Z(node2, 1)];
%             isPresent = find(ismember(points_new, midpoint, 'rows'));
% 
%             if (isempty(isPresent))
%                 points_new = [points_new; midpoint];
%                 structure_(i, ksf(j) + 1) = size(points_new, 1) + sizeOfPnts_now;
%             else
% 
%                 if (size(isPresent, 1) == 1)
%                     structure_(i, ksf(j) + 1) = isPresent(1, 1) + sizeOfPnts_now;
%                 else
%                     error('Duplicate points exist in array of points_new 2')
%                 end
% 
%             end
% 
%             clear midpoint
%         end
% 
%         % --------------------------- lateral
%         exs = 17;
% 
%         for j = 1:4
%             node1 = structure_(i, ksd(j));
%             node2 = structure_(i, ksd(j) + 8);
%             
%             midpoint = 0.5 * [X(node1, 1) + X(node2, 1), ...
%                             Y(node1, 1) + Y(node2, 1), ...
%                             Z(node1, 1) + Z(node2, 1)];
% 
%             isPresent = find(ismember(points_new, midpoint, 'rows'));
% 
%             if (isempty(isPresent))
%                 points_new = [points_new; midpoint];
%                 structure_(i, exs) = size(points_new, 1) + sizeOfPnts_now;
%             else
% 
%                 if (size(isPresent, 1) == 1)
%                     structure_(i, exs) = isPresent(1, 1) + sizeOfPnts_now;
%                 else
%                     error('Duplicate points exist in array of points_new')
%                 end
% 
%             end
% 
%             clear midpoint
%             exs = exs + 1;
%         end
% 
%     end
%     Vpoints = [Vpoints; points_new];
end
