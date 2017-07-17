function mdlsShowBlacklist(D,labels,ct,outdir)
    if ~exist('outdir','var')
        outdir = '.';
    end
    Db = [];
    K  = [];
    C = length(D);
    %
    % Really damn big dictionary
    %
    c = 1;
    for i1 = 1:C
        for i2 = 1:C
            if i1==i2
                c = c + 1;
                continue;
            end
            G       = abs(D{i2}'*D{i1});
            [mv,mi] = max(G);
            max(mv)
            blist   = find(mv >= ct);
            fprintf('Blacklist between dics %d and %d: %d\n',labels(i1),labels(i2),length(blist));
            if ~isempty(blist)
                Dbi     = D{i1}(:,blist);
                %mdlsFigure( sprintf(...
                %    'Blacklist between dics %d and %d',i1,i2) );
                mdlsFigure('Blacklist');
                subplot(C,C,c); 
                I = mdlsDictDisplay(Dbi);
                imagesc(I); 
                imwrite(I,sprintf('%s/blacklist_%d_vs_%d.png',outdir,labels(i1),labels(i2)));
                colormap gray; axis off;
                Db      = [Db Dbi];
            end
            c = c + 1;
        end
    end    
end