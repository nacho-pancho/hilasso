function s=mdlsPrintLatexTable(table,varargin)

    [nrows,ncols]=size(table);
    
    markBest='none';
    critBest='max';
    markWorst='none';
    critWorst='min';
    
    justif=repmat('c',1,ncols);
    precSpec = cell(1,ncols);
    for j = 1:ncols
        precSpec{1,j} = '%f';
    end
    %
    % Parse parameters
    %
    for i = 1:2:length(varargin)
        if isequal(varargin{i}, 'RowLabel')
            justif=['l|' justif];
            rowLabel=varargin{i+1};
        elseif isequal(varargin{i}, 'PrecSpec') % precision fields in
                                                % printf format
            precSpec=varargin{i+1};
        elseif isequal(varargin{i}, 'ColLabel')
            colLabel=varargin{i+1};
        elseif isequal(varargin{i}, 'MarkBest')
            markBest=varargin{i+1};
        elseif isequal(varargin{i}, 'critBest')
            bestCrit=varargin{i+1};
        elseif isequal(varargin{i}, 'MarkWorst')
            markBest=varargin{i+1};
        elseif isequal(varargin{i}, 'critWorst')
            bestCrit=varargin{i+1};
        else
            error(['mdlsPrintLatexTable: Unknown option : ' varargin{i}]);
        end
    end    
    %
    % Table Header
    %
    if exist('colLabel','var')
        s = [ mdlsJoinCell(colLabel,' & ') '\\ \hline' sprintf('\n')];
        if exist('rowLabel','var')
            s = [ '     & ' s ];
        end
    end
    %
    % Body
    %
    for i=1:nrows
        if exist('colLabel','var')
            srow = [colLabel{i} ' & ']; 
        else
            srow = '';
        end
        for j=1:ncols-1
            srow = [ srow sprintf(precSpec{1,j}, table(i,j)) ' & '];
        end
        fmt = [ '%s' precSpec{1,ncols} '\\\\\n' ];
        s = [ s sprintf(fmt, srow, table(i,ncols)) ];
    end
    %
    % put it in a tabular environment
    %
    s = latexBlock('tabular',s,['{' justif '}']);
    s = [ s '\caption{\label{tab:noname}}' ];
    %
    % put it in a table environment
    %
    s = latexBlock('table',s);
end

function s=latexBlock(block_name,block_content,block_args) 
  if nargin < 3
    s=sprintf('\\begin{%s}\n%s\n\\end{%s}\n',  block_name, block_content,block_name);
  else
    s=sprintf('\\begin{%s}%s\n%s\\end{%s}\n',  block_name, block_args,block_content,block_name);
  end
end