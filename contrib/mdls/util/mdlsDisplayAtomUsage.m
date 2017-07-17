% function I = mdlsDisplayAtomUsage(usageStats,maxVal,resolution,mode)
%
% Category: data display
%
% Description:
%
% Computes an image that shows the usage of each atom of a dictionary by each class
% in a set of classes in a very compact way.
% 
% Input:
% usageStats ........ (1xK)  contains the
%                     usage of each of the K atoms by a set of data samples.
%                     The usage of each atom is computed as the sum of
%                     the absolute value of the reconstruction coefficients
%                     associated to that atom.         
% maxVal ............ Maximum value, used for normalization
% resolution ........ Number of pixels in the scale, gives the height of the output matrix.
%
% Output:
% I ................. Image that shows the atom usage graphically.
%
% @author Ignacio Ramirez and Federico Lecumberry ({nacho,fefo}@fing.edu.uy)
%
function I = mdlsDisplayAtomUsage(usageStats,maxVal,resolution,mode)
  if ~exist('resolution','var')
    resolution = 300;
  end
    
  if ~exist('maxVal','var')
    maxVal = max(usageStats);
  end

  if ~exist('mode','var')
    mode='bars';
  end

  m = resolution; 
  n = length(usageStats);
  
  I=ones(m,n);
  if isequal(mode,'bars')
    m2 = ceil(m/2);
    C = m2/maxVal;
    for k=1:n
      v=usageStats(k);
      m0=m2-ceil(C*v);
      m1=m2+ceil(C*v);
      if (m0 < 1)
        m0 = 1;
      end
      if (m1 > m)
        m1 = m;
      end
      I(m0:m1,k)=0;
    end
  else
    C = 1/maxVal;
    for k=1:n
      v = C * usageStats(k);
      I(:,k) = v;
    end
  end
end
