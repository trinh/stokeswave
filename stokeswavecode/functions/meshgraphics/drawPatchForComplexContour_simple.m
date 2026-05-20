function handles = drawPatchForComplexContour_simple(Cfull, ratio_normal, reimname)
% DRAW_PATCHES draws a series of patches along a curve

% We have to split up the surface each time a point repeats
num_of_cuts = 1;

h_patch_right = cell(1, num_of_cuts);
h_patch_left = cell(1, num_of_cuts);
h_edge_right = cell(1, num_of_cuts);
h_edge_left = cell(1, num_of_cuts);
for j = 1:num_of_cuts
    
    if num_of_cuts == 1
        C = Cfull;
    else
        if j == 1
            C = Cfull(:, 1:ind(1));
        elseif j == num_of_cuts
            C = Cfull(:, ind(j-1) + 1:end);
        else
            C = Cfull(:, ind(j-1) + 1:ind(j));
        end
    end
    
    pp = getPatchForComplexContour(C, ratio_normal, reimname);
    
    
    verts_right = pp.verts_right;
    faces_right = pp.faces_right;
    verts_left = pp.verts_left;
    faces_left = pp.faces_left;
    
    h_patch_right{j} = patch('Vertices', [real(verts_right(:, 1:2)), verts_right(:,3)], 'Faces', faces_right);
    set(h_patch_right{j}, 'EdgeColor', 0.9*[1, 1, 1], 'FaceColor', 'w');
    hold on
    
    h_patch_left{j} = patch('Vertices', [real(verts_left(:, 1:2)), verts_left(:,3)], 'Faces', faces_left);
    set(h_patch_left{j}, 'EdgeColor', 0.9*[1, 1, 1], 'FaceColor', 'w');
    
    A2 = verts_right(faces_right(:,2), :);
    A3 = verts_right(faces_right(:,3), :);
    x_edge = [A2(:,1), A3(:,1)];
    y_edge = [A2(:,2), A3(:,2)];
    f_edge = [A2(:,3), A3(:,3)];
    h_edge_right{j} = plot3(x_edge, y_edge, f_edge, 'k', 'LineWidth', 2);
    
    A2 = verts_left(faces_left(:,2), :);
    A3 = verts_left(faces_left(:,3), :);
    x_edge = [A2(:,1), A3(:,1)];
    y_edge = [A2(:,2), A3(:,2)];
    f_edge = [A2(:,3), A3(:,3)];
    h_edge_left{j} = plot3(x_edge, y_edge, f_edge, 'k', 'LineWidth', 2);
    
    handles = [];
    handles.patch_right{j} = h_patch_right{j};
    handles.patch_left{j} = h_patch_left{j};
    handles.edge_right{j} = h_edge_right{j};
    handles.edge_left{j} = h_edge_left{j};
    
    hold off
    
end

