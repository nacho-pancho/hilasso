function mdlsPrintByName(figdir,figname,fig_handle_name)
  if ~exist('fig_handle_name','var')
    fig_handle_name=figname; % very reasonable default
  end
  f=gcf(); % push
  figByName(fig_handle_name);
  print('-depsc2','-r300','-tiff',[figdir '/' figname '.eps']);
  print('-dpng','-r300',[figdir '/' figname '.png']);
  figure(f); % pop
end
