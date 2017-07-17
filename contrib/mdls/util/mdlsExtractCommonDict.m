function D2=mdlsExtractCommonDict(D,maxmu)
if ~iscell(D)
  error('D must be a cell of at least 2 dictionaries');
end
NC=length(D);
if NC < 2
  error('D must have at least 2 dictionaries');
end

if ~exist('maxmu','var')
maxmu = 0.95;
end

if ~exist('maxclashes','var')
maxclashes = 0.9; % 90% of other atoms
end

%
% scan dictionaries
%
K = zeros(1,NC);
Kt = zeros(1,NC+1);
Kt(1) = 1;
for g=1:NC
    K(g) = size(D{g},2);
    Kt(g+1) = Kt(g) + K(g);
end
%
% concatenated dictionary
%
Dt = zeros(size(D{1},1),Kt(end)-1);
for g=1:NC
    Dt(:,Kt(g):(Kt(g+1)-1)) = D{g};
end
%
% S counts, for each atom in Dt, how many colinear
% atoms there are in each other dictionary.
%
S=zeros(size(Dt,2),NC);
for g=1:NC
  S(:,g) = sum(abs(Dt'*D{g}) > maxmu, 2);
end
St = abs(Dt'*Dt) > maxmu;
%
% find those atoms for which there is another very
% similar one in each other dictionary
%
common_atoms = find(sum(S>0,2) == NC); 
%
% create common dictionary
%
D2 = cell(1,NC+1);
D2{NC+1} = Dt(:,common_atoms);
% 
% remove the common atoms, as well as those that
% are very close to them, from the dictionaries
fprintf('Found %d common atoms\n',length(common_atoms));
for c=1:length(common_atoms)
  common_atoms = union(common_atoms,find(St(common_atoms(c),:)));
end
fprintf('Pruning %d atoms in total\n',length(common_atoms));
%
% remove common atoms from original dictionaries
%
for g=1:NC
    idx_common_atoms_in_g = (Kt(g) <= common_atoms) & (common_atoms < Kt(g+1));
    common_atoms_in_g = common_atoms(idx_common_atoms_in_g) - (Kt(g) + 1);
    D2{g} = D{g}( :, setdiff(1:K(g),common_atoms_in_g) );                      
end
end