% function Y = col_group_vector_soft(X,tau,groups)
% 
% INPUT
% X ............. (KxN)
% groups ........ (Kx1) vector indicating the group stucture
%                 for the group-vector-soft thresholding operation
%                 Each set of elements that are mutually equal 
%                 define a group. 
%                 Example: if groups = [1 1 1 2 2 2 3 3 2 2 1 1 3 3],
%                 there are 3 groups and 
%                   X([1,2,3,11,12],:) are in group 1,
%                   X([4,5,6,9,10],:) are in group 2,
%                   X([7,8,13,14],:) are in group 3.
%                 All elements of groups must be integers in the range 
%                  1...num_groups
% tau ........... threshold penalty
%
% OUTPUT
% Y ............. thresholded X
%
function Y = col_group_vector_soft(X,tau,groups)

  num_groups = max(groups);
  if min(groups)~=1
     error(['Wrong group structure vector'])
  end
  Y = zeros(size(X));
  for i=1:num_groups
      thisgroup = find(groups == i);
      Y(thisgroup,:) = matrix_soft(X(thisgroup,:),tau);
  end
end
