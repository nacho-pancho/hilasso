function mdlsPrint(figname,figdir)
  if ~exist('figdir','var')
    figdir='.';
  end
  print('-depsc2','-r300','-tiff',[figdir '/' figname '.eps']);
  print('-dpng','-r300',[figdir '/' figname '.png']);
end
