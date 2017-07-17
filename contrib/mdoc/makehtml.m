function [am,lab,longlab,dirs]=makehtml(list,ttl,tmpl,itmpl,dir,ignore,require)
%function [am,lab,longlab,dirs]=makehtml(list,ttl,tmpl,itmpl,dir,ignore,require)
%
% $Id: makehtml.m,v 1.2 2004-03-25 18:33:49-06 brinkman Exp brinkman $
% (c) 2002 Peter Brinkmann, brinkman@math.uiuc.edu
%
% Usage example: makehtml('makehtml','makehtml utility');
%
% Parameters:
%	list: string or cell array of strings containing list of m-files
%		to be documented
%	ttl: (optional) title; defaults to ''
%	tmpl: (optional) template file for module pages,
%		defaults to 'mdoctemplate.html'
%	itmpl: (optional) template file for index page,
%		defaults to tmpl
%	dir: (optional) directory where genhtml creates its output;
%		defaults to 'content/'
%	ignore, require: (optional) see layout.m
%
% Returns:
%	see depgraph.m

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

	disp('$Id: makehtml.m,v 1.2 2004-03-25 18:33:49-06 brinkman Exp brinkman $')
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
		tmpl='mdoctemplate.html';
	end;
	if nargin<4
		itmpl=tmpl;
	end;
	if nargin<5
		dir='content/';
	end
	if nargin<6
		ignore={};
	end
	if nargin<7
		require={pwd};
	end

	[am,lab,longlab,dirs]=depgraph(list,ignore,require);
	genhtml(am,lab,ttl,tmpl,itmpl,dir);

