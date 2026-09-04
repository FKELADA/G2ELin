function [rightV,leftV] = EigVec(A)
    [V,D,W] = eig(double(vpa(subs(A)))); % obtain eigenvalues and eigenvectors

    % Eigenvectors
    % columns of V are the right eigenvectors
    % columns of W are the left eigenvectors
    coef = diag(W'*V); % for normalization of eigenvectors
    rightV = V;
    leftV = zeros(size(A));
    for i = 1:length(A)
        leftV(i,:) = 1/coef(i)*(W(:,i)'); % normalization of left eigenvectors
    end
    % columns of rightV are the right eigenvectors
    % rows of leftV are the left eigenvectors
    % rightV and leftV are orthonormalized -> leftV*rightV = I
    % ATTENTION : leftV is not like Matlab W, it is actually W' normalized;
    % the ~rows~ of leftV are eigenvectors, not the columns
end