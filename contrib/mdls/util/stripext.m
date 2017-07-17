function fname=stripext(object_id)
  extindex = rindex(object_id,".");
  if extindex > 1
    fname = substr(object_id,1,extindex-1);
  else
    fname= object_id;
  end
endfunction
