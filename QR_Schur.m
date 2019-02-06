% Algorithm to reduce the matrix A to its Schur form
%% variable declarations
A = randi(3,3)

[Q, H] = hess(A);
localH = H;
iterations = 0;
m =  length(A);
n = m
tolerance = 1.5;

%% main loop
while m>1
    iterations
    Q
    localH
    %% get the element to reduce
    target = localH(m,m);
    %% compute Qm and R
    [R, Qm] = given (localH(1:m,1:m) - target .* eye(m));
    %% update localH
    localH(1:m,1:m) =  R * Qm + target .* eye(m)
    
    if m<n
        localH(1:m, m+1:n) = Qm' * localH(1:m, m+1:n);
    end

    %% update Q
    Q(:,1:m) = Q(:,1:m) * Qm;
    %% mark iteration and check if tolerance is met
    iterations = iterations+1;
    if abs(localH(m,m-1)) < tolerance
        localH(m,m-1) = 0;
        m = m-1;
    end
end

%% check
Q * localH * Q'