

% D{1} = getDict('data/brodatz/D1.png',50,10,9,300);
% D{2} = getDict('data/brodatz/D52.png',50,10,9,300);
% D{3} = getDict('data/brodatz/D111.png',50,10,9,300);
% D{4} = getDict('data/brodatz/D49.png',50,10,9,300);
% 
% for i=1:4
%     Do{i} = D{i}{1};
% end

% [Y{1},I{1}] = getTestImage('data/brodatz/D1.png',50,80,10,9);
% [Y{2},I{2}] = getTestImage('data/brodatz/D52.png',50,80,10,9);
% [Y{3},I{3}] = getTestImage('data/brodatz/D111.png',50,80,10,9);
% [Y{4},I{4}] = getTestImage('data/brodatz/D49.png',50,80,10,9);


lambda1 = [170 200 250];

lambda2 = [1 3 5];

for i=1:length(lambda1)
    
    for j=1:length(lambda2)
        
%         for r=1:4
                disp([i j r])
%             if sum([i,j,r]==[1,1,1])~=3 && sum([i,j,r]==[1,1,2])~=3 && sum([i,j,r]==[1,1,3])~=3
        [Xr2{i}{j},A2{i}{j}] = HIlassoColMethod(Yr,Do,lambda1(i),lambda2(j),0.01);
%             end
        e = 0;
        for h=1:4
        Io = mdlsReconstruct(Xr2{i}{j}{h},size(Io,1),size(Io,2),9);
        e = e+norm(Io-I{h},'fro')^2/norm(Ie,'fro')^2;
        end        
        error(i,j) = e;
        error
        aset(i,j) = mean(sum(A2{i}{j}~=0));
        pause(0.5)

        
        
    end
end
        

% lambdaL = 50:10:180;
% 
% 
% for i=1:length(lambdaL)
% 
% 
%     disp(lambdaL(i))
%     [XrL{i},AL{i}] = lassoMethod(Yr,Do,lambdaL(i));
%     e = 0;
%     for h=1:4
%         Io = mdlsReconstruct(XrL{i}{h},size(Io,1),size(Io,2),9);
%         e = e+norm(Io-I{h},'fro')^2/norm(Ie,'fro')^2;
%     end
%     errorL(i,j) = e;
%     errorL
%     asetL(i,j) = mean(sum(AL{i}~=0));
%     pause(0.5)
% 
% 
% end
