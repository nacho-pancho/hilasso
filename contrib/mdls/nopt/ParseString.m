function out = ParseString(string,delimiter,key)
  % out = ParseString(string,delimiter,key)

  if ~isempty(strfind(string,delimiter))
    rem = string;
    k = 1;
    while ~isempty(rem)
      [val,rem] = strtok(rem,delimiter);
      if ~isempty(val)
        out(k) = {val};
        k = k + 1;
      end
    end

    if isequal(nargin,3) && strcmpi(key,'last')
      out = out{end};
    end
    if isequal(nargin,3) && strcmpi(key,'first')
      out = out{1};
    end
  else
    out = string;
  end
end
