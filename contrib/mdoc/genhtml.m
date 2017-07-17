function genhtml(am,lab,ttl,tmpl,indtmpl,dir)
% function genhtml(am,lab,ttl,tmpl,indtmpl,dir)
%
% $Id: genhtml.m,v 1.17 2003-01-05 05:37:05-06 brinkman Exp $
% (c) 2002 Peter Brinkmann, brinkman@math.uiuc.edu
%
% genhtml: generates html documentation from a dependency graph and labels;
%	see also drawlayout.m
%
% Usage example: genhtml([1,1;0,1],{'module1','module2'})
%
% If N is the number of vertices in the dependency graph, then genhtml
% creates N+1 content files for Genpage; one index file and one file
% for each vertex.
%
% Each content file has four content tags; 'title' for the title of
% the page, 'helptext' containing the result of help(module), and
% 'dependencies' listing links to html-files for modules that the
% current modules depends on, and 'invdependencies' listing links
% to modules that depend on the current module.
%
% The dependencies are listed as rows of tables with one cell per row.
% The Genpage template file needs to enclose the list of dependencies
% in table tags.
%
% You can also choose separate layout files for the index file and the
% files for individual modules.
%
% Parameters:
%	am,lab,ttl: adjacency matrix, names of modules, optional title
%		(see drawlayout.m for details)
%	dir: (optional) html directory, defaults to './content/'
%	tmp: (optional) genpage template file, defaults to 'mdoclayout.html'
%	indtmp: (optional) genpage template file for index page, defaults to tmpl

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

	if nargin<3
		ttl='';
	end
	if nargin<4
		tmpl='mdoclayout.html';
	end
	if nargin<5
		indtmpl=tmpl;
	end
	if nargin<6
		dir='content/';
	end

	N=length(lab);

	fid=myfopen(dir,'index.content');

	fprintf(fid,['%%template%%\n' indtmpl '\n%%title%%\n' ttl ...
		'\n%%helptext%%\n\n%%dependencies%%\n']);
	
	for i=1:N
		fprintf(fid,'<tr><td>%s</td></tr>\n',link(lab{i}));
	end

	fclose(fid);

	[ii,jj]=meshgrid(1:N);
	depcnt=sum((ii~=jj) & am,2);
	invdepcnt=sum((ii~=jj) & am,1);

	for i=1:N
		fid=myfopen(dir,[filename(lab{i}) '.content']);
		fprintf(fid,['%%template%%\n' tmpl '\n%%title%%\n' lab{i}...
			'\n%%helptext%%\n' help(lab{i}) '\n%%dependencies%%\n']);
		if depcnt(i)
			for j=1:N
				if i~=j & am(i,j)
					fprintf(fid,'<tr><td>%s</td></tr>\n',link(lab{j}));
				end
			end
		else
			fprintf(fid,'<tr><td>---</td></tr>\n');
		end
		fprintf(fid,'%%invdependencies%%\n');
		if invdepcnt(i)
			for j=1:N
				if i~=j & am(j,i)
					fprintf(fid,'<tr><td>%s</td></tr>\n',link(lab{j}));
				end
			end
		else
			fprintf(fid,'<tr><td>---</td></tr>\n',link(lab{j}));
		end
		fclose(fid);
	end


function s=link(fn)
	s=['<a href="' filename(fn) '.html">' fn '</a>'];


function s=filename(fn)
	s=['mdoc_' fn];		% need to prefix the module name with something
		% in case there's a module called 'index.m'


function f=myfopen(dir,fn)
	[f,msg]=fopen([dir '/' fn],'w');
	if f<0
		error(msg);
	end

