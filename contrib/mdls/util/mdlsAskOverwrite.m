function reply = AskOverwriteFile(filename,backup)
if exist(filename,'file')
  if exist('backup','var') && backup
    comando = sprintf('mv %s %s.prev',filename,filename);
    system(comando);
    reply = true;
  else
    texto = sprintf('%s exists. Overwrite it? [y|N] ',filename);
    reply = input(texto,'s');
    if isequal(reply,'Y') || isequal(reply,'y')
      reply = true;
      fprintf(1,'\b (yes)\n')
    else
      reply = false;
      fprintf(1,'\b (no)\n')
    end
  end
else
  reply = true;
end
