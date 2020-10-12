function [pval, tastepairs] = fun_compare_tastes_v2(n,comparison, binnedC, tastes)
tastepairs=[];
switch n
    
    case 1
        x=length(tastes)/2;
        pval = NaN( size(binnedC,4),x);
        for j=1:x
            for i=1:size(binnedC,4)
                
                if j ==1
                    pval(i,j)=ranksum(binnedC(:,comparison,1,i),binnedC(:,comparison,length(tastes),i));
                end
                if j ==2
                    pval(i,j)=ranksum(binnedC(:,comparison,2,i),binnedC(:,comparison,length(tastes)-1,i));
                end
                if j ==3
                    pval(i,j)=ranksum(binnedC(:,comparison,3,i),binnedC(:,comparison,length(tastes)-2,i));
                end
            end
        end
        
        
    case 2
        x=4;
        pval = NaN( size(binnedC,4),x);
        for j=1:x
            for i=1:size(binnedC,4)
                
                if j ==1
                    pval(i,j)=ranksum(binnedC(:,comparison,1,i),binnedC(:,comparison,3,i));
                end
                if j ==2
                    pval(i,j)=ranksum(binnedC(:,comparison,2,i),binnedC(:,comparison,4,i));
                end
                
            end
        end
    case 3
        
        pval = NaN( size(binnedC,4),1);
        
        for i=1:size(binnedC,4)          
            pval(i,1)=ranksum(binnedC(:,comparison,1,i),binnedC(:,comparison,2,i));
        end
        
    case 4
        tastearray=1:1:length(tastes);
        tastepairs=nchoosek(tastearray,2);
        pval = NaN(size(binnedC,4),size(tastepairs,1));
%         if length(tastepairs)>1
        for j=1:size(tastepairs,1)
            for i=1:size(binnedC,4)
            
                    pval(i,j)=ranksum(binnedC(:,comparison,tastepairs(j,1),i),binnedC(:,comparison,tastepairs(j,2),i));
                      
            end
        end
        
end
end