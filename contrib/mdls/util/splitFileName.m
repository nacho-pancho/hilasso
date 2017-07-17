function [fpath,fname,fext]=splitFileName(full_path)
  dot_pos = rindex(full_path,".");
  if dot_pos == 0
    fext = "";
    fname=full_path;
  else
    fext=substr(full_path,dot_pos+1);
    fname=substr(full_path,1,dot_pos-1);
  end
  path_pos = rindex(fname,"/");
  if path_pos == 0
    fpath = "./";
  else
    fpath = substr(fname,1,path_pos);
    fname = substr(fname,path_pos+1);
  end
endfunction