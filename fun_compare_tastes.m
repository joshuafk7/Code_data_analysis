function pval = fun_compare_tastes(n,comparison, binnedC, neurons, tastes)
switch n
    case 1
        x=length(tastes)/2;
        pval = NaN( length(neurons),x);
        for j=1:x
            for i=1:length(neurons)
                
                if j ==1
                    pval(i,j)=ranksum(binnedC(:,comparison,1,i),binnedC(:,comparison,6,i));
                end
                if j ==2
                    pval(i,j)=ranksum(binnedC(:,comparison,2,i),binnedC(:,comparison,5,i));
                end
                if j ==3
                    pval(i,j)=ranksum(binnedC(:,comparison,3,i),binnedC(:,comparison,4,i));
                end
            end
        end
        
        
    case 2
        x=4;
        pval = NaN( length(neurons),x);
        for j=1:x
            for i=1:length(neurons)
                
                if j ==1
                    pval(i,j)=ranksum(binnedC(:,comparison,1,i),binnedC(:,comparison,2,i));
                end
                if j ==2
                    pval(i,j)=ranksum(binnedC(:,comparison,3,i),binnedC(:,comparison,4,i));
                end
                
            end
        end
end
end
