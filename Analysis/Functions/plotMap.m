function [out] = plotMap(geom_frame, map)
%PLOTMAP plotting a 2D map
figure;
imagesc(geom_frame.xEdges, geom_frame.yEdges, map);
axis xy;                 % Keep Y axis increasing upwards
xlabel('X position (mm)');
ylabel('Y position (mm)');
colorbar;
colormap(jet); 

end

