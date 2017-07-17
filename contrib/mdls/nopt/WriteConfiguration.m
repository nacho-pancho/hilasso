function WriteConfiguration(options)
fn = fieldnames(options);
LogFile(options,'write','{\n');
for m = 1:length(options)
  for n = 1:length(fn)
    fN = fn{n};
    tmp = options(m).(fN);
    if (isnumeric(tmp) || islogical(tmp))
      estrin = num2str(tmp);
    elseif isstruct(tmp)
      options(m).(fN)(1).logFile = options.logFile;
      WriteConfiguration(options.(fN))
    else
      estrin = tmp;
    end
    LogFile(options,'write',sprintf('\t# %s: %s\n',fN,estrin));
  end
end
LogFile(options,'write','}\n');
