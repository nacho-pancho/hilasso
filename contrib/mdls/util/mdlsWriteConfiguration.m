function mdlsWriteConfiguration(options)
fn = fieldnames(options);
fprintf('{\n');
for m = 1:length(options)
  for n = 1:length(fn)
    fN = fn{n};
    tmp = options(m).(fN);
    if (isnumeric(tmp) || islogical(tmp))
      estrin = num2str(tmp);
    elseif isstruct(tmp)
      mdlsWriteConfiguration(options.(fN));
    elseif iscell(tmp)
        % skip cells for now
        estrin = sprintf('(%dx%d cell)',size(tmp,1),size(tmp,2));
    else
      estrin = tmp;
    end
    fprintf('\t# %s: %s\n',fN,estrin);
  end
end
fprintf('}\n');
