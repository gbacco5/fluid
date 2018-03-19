function mi_draw_arc(x1, y1, x0, y0, x2, y2, maxseg)

R1 = hypot(x1-x0, y1-y0);
R2 = hypot(x2-x0, y2-y0);

if round(R1*1e6)/1e6 ~= round(R2*1e6)/1e6
  error('Incompatible coordinates for an arc.')
else
  R = R1;
end

alpha2 = atan2(y2-y0,x2-x0);
alpha1 = atan2(y1-y0,x1-x0);
alpha = alpha2 - alpha1;

mi_drawarc(x1,y1, x2,y2, alpha*180/pi, maxseg)
  
end