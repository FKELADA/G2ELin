function [partMatrix] = PartMatrix(A,rightV,leftV,ylabels,xlabels,n_exc_ln_ld,start_from_eig,end_eig)
    partMatrix = zeros(size(A));
    den = zeros(length(A));
    for j = 1:length(A)
        sum = 0;
        for k = 1:length(A)
            sum = sum + abs(leftV(j,k))*abs(rightV(k,j));
        end
        den(j) = sum;
    end
    for i = 1:length(A)
        for j = 1:length(A)
            partMatrix(i,j) = abs(rightV(i,j))*abs(leftV(j,i))/den(j);
            % rightV(i,j) is the activity of ith state in jth mode
            % leftV(j,i) is the contribution of this activity to jth mode
            % P(i,j) is the relative participation of ith state in jth mode
            % each row of P corresponds to a state
            % each column of P corresponds to a mode
        end
    end
    
    figure
    heatmap(xlabels,ylabels,partMatrix.*100)

end