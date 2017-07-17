function [am,lab,longlab,dirs,fig]=makefig(list,ttl,fn,ignore,require)
% function [am,lab,longlab,dirs]=makefig(list,ttl,fn,ignore,require)
%
% $Id: makefig.m,v 1.4 2003-05-30 04:15:07-05 brinkman Exp $
% (c) 2002 Peter Brinkmann, brinkman@math.uiuc.edu
%
% makefig: tool for creating a documentation GUI for Matlab packages
%
% Usage example: makefig('makefig','makefig utility','foo.fig')
%
% Parameters:
%	list: string or cell array of strings containing list of m-files
%		to be documented
%	ttl: (optional) title; defaults to ''
%	fn: (optional) file name for saving figure; defaults to []
%	ignore, require: (optional) see graphlayout.m
%
% Returns:
%	see depgraph.m, plus fig, which is the handle of the figure

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

	disp('$Id: makefig.m,v 1.4 2003-05-30 04:15:07-05 brinkman Exp $')
	disp('Copyright (c) 2001, Peter Brinkmann (brinkman@math.uiuc.edu)')
	disp('This program is distributed in the hope that it will be useful,')
	disp('but WITHOUT ANY WARRANTY; without even the implied warranty of')
	disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the')
	disp('GNU General Public License for more details.')
	disp('');

	if nargin<2
		ttl='';
	end
	if nargin<3
		fn=[];
	end
	if nargin<4
		ignore={};
	end
	if nargin<5
		require={pwd};
	end

	[am,lab,longlab,dirs]=depgraph(list,ignore,require);
	xyc=graphlayout(am);
	fig=drawlayout(xyc,am,lab,ttl);

	if fn
		saveas(gcf,fn);
	end

