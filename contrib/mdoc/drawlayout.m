function fig=drawlayout(xyc,am,lab,ttl)
% function fig=drawlayout(xyc,am,lab,ttl)
%
% $Id: drawlayout.m,v 1.13 2003-05-30 04:20:17-05 brinkman Exp $
% (c) 2002 Peter Brinkmann, brinkman@math.uiuc.edu
%
% drawlayout: draws a graph given the coordinates and labels of vertices
%
% Usage example: drawlayout([1,2;3,1],[0,1;0,0],{'a','b'},'title');
%
% Parameters
%	xyc: 2xN matrix; the first row is the list of x-coordinates, and the
%		second row is the list of y-coordinates
%	am: NxN adjacency matrix of graph
%	lab: (optional) cell array of N labels
%	ttl: (optional) title of figure
%
% Returns
%	fig: handle of graphics window

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

	fig=figure;

	cla;
	axis manual;
	axis([0,1,0,1]);
	hold on;
	set(gca,'xtick',[]);
	set(gca,'ytick',[]);

	N=length(am);
	xx=xyc(1,:);
	yy=xyc(2,:);

	xx=xx-min(xx);
	yy=yy-min(yy);
	xx=0.1+0.8*xx/max(1,max(xx));
	yy=0.1+0.8*yy/max(1,max(yy));

	boxes={};
	if nargin>2
		for i=1:N
			h=text(xx(i),yy(i),lab{i},...
				'Interpreter','none',...
				'HorizontalAlignment','center',...
				'ButtonDownFcn','doctool(gcbo)');
			boxes{i}=get(h,'Extent');
		end
	else
		for i=1:N
			plot(xx(i),yy(i),'o');
			boxes{i}=[xx(i),yy(i),0,0];
		end
	end

	for i=1:N
		for j=1:N
			if i~=j & am(i,j)
				connect(boxes{i},boxes{j});
			end
		end
	end

	if nargin>3
		title(ttl);
	end


function connect(b0,b1)
	p0=[b0(1)+b0(3)/2,b0(2)+b0(4)/2];
	p1=[b1(1)+b1(3)/2,b1(2)+b1(4)/2];
	r=p1-p0;

	t0=min(abs([b0(3)/2,b0(4)/2]./r));
	t1=min(abs([b1(3)/2,b1(4)/2]./r));
	p0=p0+t0*r;
	p1=p1-t1*r;

	r=r/sqrt(r(1)*r(1)+r(2)*r(2));
	phi=pi/8;
	al=p1-0.01*r*[cos(phi) sin(phi);-sin(phi) cos(phi)];
	ar=p1-0.01*r*[cos(phi) -sin(phi);sin(phi) cos(phi)];

	line([p0(1),p1(1)],[p0(2),p1(2)]);
	line([p1(1),al(1)],[p1(2),al(2)]);
	line([p1(1),ar(1)],[p1(2),ar(2)]);

