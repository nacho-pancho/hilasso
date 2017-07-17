%
% if some atoms in a dictionary D are found to be unused (according to some threshold
% relative to the average coefficient value.
% they are replaced by a noisy average of the L most used atoms
%
function D = ReviveAtoms(D,A,thres,L)
  if ~exist('thres','var')
    thres = 1e-4;
  end
  if ~exist('L','var')
   L = 5;
  end
   K =  size(D,2);
   A=abs(A);
   S=sum(A,2);
   ref = max(mean(A,2)) * thres;
   [SS,Si]=sort(S);
   deadAtomInd  = Si ( SS < ref );
   aliveAtomInd = Si (SS >= ref );
   sigma = sqrt(max(sum(D.*D))) * 0.05;
   layers = ceil(length(aliveAtomInd) * 0.1);
   if ~isempty(deadAtomInd)
     D(:,deadAtomInd) = GeneratePatchesDictionary(D,length(deadAtomInd),sigma,layers);
   end  
end
