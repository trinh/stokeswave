function handles = drawPatchForComplexContour2(Cfull, ratio_normal, reimname)
% DRAW_PATCHES draws a series of patches along a curve

% We have to split up the surface each time a point repeats
Cdiff_xy = diff(Cfull(1,:));
ind = find(Cdiff_xy==0);

% Remove duplicates
Cfull(:, ind+1) = [];
C = Cfull;

pp = getPatchForComplexContour(C, ratio_normal, reimname);

% We now re-assemble, but remove the 

j = 1;
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
keyboard

end