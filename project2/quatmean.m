function q_mean = quatmean(Y)
    % taking the principle eigenvector of the mean fo the outter product of
    % the quaternions. (explained in class)
    M = (1/(2*n))*(Y*Y'); %[4x4]
    [V,D] = eig(M);
    [max_lambda, ind] = max(diag(D));
    q_mean = V(:,ind);

end

