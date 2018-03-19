function draw_fluid_barrier(b)

for bkk = 1:length(b)
  xE = b(bkk).X(1);
  yE = b(bkk).Y(1);
  xEOC = b(bkk).X(2);
  yEOC = b(bkk).Y(2);
  xC = b(bkk).X(3);
  yC = b(bkk).Y(3);
  
  xD = b(bkk).X(end-2);
  yD = b(bkk).Y(end-2);
  xDOE = b(bkk).X(end-1);
  yDOE = b(bkk).Y(end-1);
  
  X = b(bkk).X(3:end-2);
  Y = b(bkk).Y(3:end-2);
  
  mi_drawpolyline([X, Y])
  
  mi_draw_arc(xE,yE, xEOC,yEOC, xC,yC, 1)
  mi_draw_arc(xD,yD, xDOE,yDOE, xE,yE, 1)
  
end

end