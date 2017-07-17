function [am,lab,longlab,dirs]=depgraph(list,ignore,require)
% function [am,lab,longlab,dirs]=depgraph(list,ignore,require)
%
% $Id: depgraph.m,v 1.8.1.13 2004-03-25 18:17:35-06 brinkman Exp $
% (c) 2002 Peter Brinkmann, brinkman@math.uiuc.edu
%
% depgraph: computes a dependency graph for a list of m-files
%
% Usage example: [am,lab]=depgraph('iode',{'/usr/local/','swapcursor'});
%
% Parameters:
%	list: m-file or cell array of m-files
%	ignore: (optional) cell array of strings; all files whose complete
%		path contains an element of ignore will be ignored; defaults to {}
%	require: (optional) cell array of strings; all files whose complete
%		path does _not_ contain any element of require will be ignored;
%		defaults to {pwd}
%
% Note that files that match the 'ignore' criterion as well as the 'require'
% criterion will be ignored. If require is empty (i.e., {}), then only paths
% matching the 'ignore' criterion will be ignored.
% 
% Returns:
%	am: adjacency matrix of dependency graph; am(i,j) is true iff
%		the i-th element of lab depends on the j-th element of lab
%	lab: short labels (base names) of files
%	longlab: long labels (full paths) of files
%	dirs: list of directories containing the files in lab

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

	if ischar(list)
		list={list};
	end
	if nargin<3
		require={pwd};
	end
	if nargin<2
		ignore={};
	end

	longlab={};
	for i=1:length(list)
		longlab{i}=which(list{i});
		if length(longlab{i})==0
			error([list{i} ' not found!']);
		end
	end
	longlab=unique(longlab);

	am=[];
	i=1;
	while i<=length(longlab)
		raw=depfun(longlab{i},'-toponly');
		longlab=append(longlab,raw,ignore, require);
		for j=1:length(raw)
			ind=find(ismember(longlab,raw{j}));
			if ind
				am(i,ind)=1;
			end
		end
		i=i+1;
	end

	lab={};
	dirs={};
	for i=1:length(longlab)
		lab{i}=basename(longlab{i});
		d=dirname(longlab{i});
		if ~ismember(d,dirs)
			dirs{length(dirs)+1}=d;
		end
	end

	dirs=sort(dirs);
	[lab,ind]=sort(lab);
	longlab=longlab(ind);
	am=am(ind,:);
	am=am(:,ind);


function ind=matches(s,lst)
	for i=1:length(lst)
		if (length(lst{i})<=length(s)) & findstr(lst{i},s)
			ind=i;
			return
		end
	end
	ind=0;


function [lm,appflag]=append(l1,l2,ignore,require)
	appflag=0;
	lm=l1;
	for i=1:length(l2)
		if ~ismember(l2{i},lm)
			if (((length(require)==0) | (matches(l2{i},require))) ...
					& (~matches(l2{i},ignore)))
				lm{length(lm)+1}=l2{i};
				appflag=1;
			end
		end
	end


function b=basename(s)
	b=s(lastslash(s)+1:end);
	m=max(findstr('.',b));
	if m
		b=b(1:m-1);
	end


function d=dirname(s)
	d=s(1:lastslash(s));


function m=lastslash(s)
% Thanks to Rajiv Narayan at Boston University for suggesting improvements!
	lst=findstr('/',s);
	if (isempty(lst))
		lst=findstr('\',s);
	end
	if (isempty(lst))
		m=0;
	else
		m=max(lst);
	end

