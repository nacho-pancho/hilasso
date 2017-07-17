function xyc=graphlayout(am,qq,km)
% function xyc=graphlayout(am,qq,km)
%
% $Id: graphlayout.m,v 1.1 2003-05-30 04:15:07-05 brinkman Exp $
% (c) 2002 Peter Brinkmann, brinkman@math.uiuc.edu
%
% graphlayout: computes a plane graphical representation of a finite graph
%
% Usage example: xyc=graphlayout(ones(5));
%
% The idea is to place charges at the vertices and strings at the edges,
% hoping to find an equilibrium position that one expects to be a good
% layout. The initial position is random, and if you're unhappy with
% the result of graphlayout, you may try running it again, hoping it'll
% arrive at a better equilibrium position this time.
%
% Parameters:
%	am: NxN adjacency matrix of graph
%	qq: (optional) vector of charges of vertices (defaults to the vector
%		of valences)
%	km: (optional) matrix of spring constants (defaults to ones(N,N))
%
% Returns:
%	xyc: 2xN matrix of coordinates [x-coordinates; y-coordinates]

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

	warning off;

	N=length(am);
	[ii,jj]=meshgrid(1:N);

	am(find(am))=1;
	am=am | am';

	if nargin<2
		qq=sum(am);
	end
	if nargin<3
		km=ones(N);
	end

	km=0.5*(km+km').*am;

	rand('seed',sum(100*clock));
	xx=N*rand(1,N);
	yy=N*rand(1,N);

	if N==0
		return
	end

	sf=1;
	i=0;
	while sf>10^-6 & i<10000
		rx=xx(ii)-xx(jj);
		ry=yy(ii)-yy(jj);
		dd=sqrt(rx.*rx+ry.*ry);	% dd(i,j) is the distance between
								% (xx(i),yy(i)) and (xx(j),yy(j))
		dd(find(dd==0))=1;		% ... avoid division by 0
	
		esfx=rx.*qq(ii).*qq(jj)./(dd.^3);	% compute electrostatic force
		esfy=ry.*qq(ii).*qq(jj)./(dd.^3);

		spfx=-rx.*km;	% compute force of strings
		spfy=-ry.*km;
	
		ffx=sum(esfx+spfx);	% total force
		ffy=sum(esfy+spfy);
	
		sf=sum(ffx.*ffx+ffy.*ffy);	% sum of the squares of lengths of forces

		xx=xx+ffx/N/5;
		yy=yy+ffy/N/5;

		i=i+1;
	end

	xyc=[xx;yy];

