function bool = in_convex_polyhedron(convex_hull, points, b, NUMthreads)
    %本程序用于判断points是否在convex_hull形成的三维凸包内
    %convex_hull为凸包或散点集
    %points为要做判断的点（nx3)
    %本程序返回类型为布尔型
    %使用实例

    %example
    %A=[0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1];
    %p=[1 2 3];
    %b=in_convex_polyhedron(A,p)

    bool = zeros(size(points, 1), 1);
    ori_set = convex_hull;
    ori_edge_index = convhull(ori_set, 'Simplify', true);
    ori_edge_index = sort(unique(ori_edge_index));

    parfor (i = 1:size(points, 1), NUMthreads)

        if (b(i) == 1) % means that this point is belonging to another half sphere
            continue
        end

        new_set = [convex_hull; points(i, :)];
        new_edge_index = convhull(new_set, 'Simplify', true);
        new_edge_index = sort(unique(new_edge_index));
        bool(i, 1) = isequal(ori_edge_index, new_edge_index);
    end

    bool = boolean(bool);
end
