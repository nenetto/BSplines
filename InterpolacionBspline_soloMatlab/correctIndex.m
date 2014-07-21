function index_new = correctIndex(index,L)

index_new = index;
indexes2 = [1:L,L-1:-1:1];
m = length(indexes2)-1;

for k = 1:length(index)
    if(L == 1)
        index_new(k) = 1;
    else
        if(index_new(k)<0)
            index_new(k) = index_new(k)+1;
        end
        i = mod(index_new(k),m);
        if(i==0)
            i=2;
        end
        index_new(k) = indexes2(i);

    end
    
end

end