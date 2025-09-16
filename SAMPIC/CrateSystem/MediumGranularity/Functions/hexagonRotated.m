function hexagonRotated(cote,x0,y0,color)
   %cote= side size;,(x0,y0) hexagon center coordinates;
   x=cote*[-0.5 0.5 1 0.5 -0.5 -1]+x0;
   y=cote*[0.866 0.866 0 -0.866 -0.866 0]+y0;
   plot(x,y,'k','Linewidth',2);grid;
   fill(x,y,color)
   axis([ x0-cote x0+cote y0-cote y0+cote]);
end