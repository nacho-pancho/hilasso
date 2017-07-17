

function plotResult(error,aSet,V,d)


for i=1:length(d) 
    if d(i) == -1
        id{i} = 1:size(error,i);
        if length(d) == 3 && i>=2
            x = i+1;
        else
            x = i;
        end
        
    elseif d(i) == 0
        if length(d) == 3 && i>=2
            param = i+1;
        else
            param = i;
        end
    else
        id{i} = d(i);
    end
end
        

color = {'b.-','r.-','k.-'};

figure
for i=1:length(V{param})
    id{param} = i;
    if i==1
        c = color{1};
    elseif i==length(V{param})
        c = color{3};
    else
        c = color{2};
    end

    subplot(1,2,1)
    hold on
    if length(d)==4
        plot(V{x},reshape(error(id{1},id{2},id{3},id{4}),1,length(V{x})),c)
    else
        plot(V{x},reshape(error(id{1},id{2},id{3}),1,length(V{x})),c)
    end
        
    
    subplot(1,2,2)
    hold on
    if length(d)==4
        plot(V{x},reshape(aSet(id{1},id{2},id{3},id{4}),1,length(V{x})),c)
    else
        plot(V{x},reshape(aSet(id{1},id{2},id{3}),1,length(V{x})),c)
    end

end