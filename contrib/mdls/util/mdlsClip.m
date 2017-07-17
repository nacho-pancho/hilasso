function y=mdlsClip(x,xmin,xmax)
y=x;
y(y>xmax)=xmax;
y(y<xmin)=xmin;
end
