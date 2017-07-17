function [f,orig_proj,stego_proj]=fisherPlotData(fld,X,Y);
  orig_proj =fld'*X;
  stego_proj=fld'*Y;
  if (mean(orig_proj)<0)
    orig_proj=-orig_proj;
    stego_proj=-stego_proj;
  end
  f=figure();
  plot(1:length(orig_proj),orig_proj,'ob;cover;',\
       1:length(stego_proj),stego_proj,'+m;stego;');
  xlabel('image');
  ylabel('FLD(moments(image))');
  title('Projection of data onto FLD vector');
end