function [Kmat_Fun,numTaus,defC0] = ParSeqModels(numC,modelType)

% Build the rate matrix functions
switch modelType
    case 'Parallel'
        Kfun_str = 'diag(-k)';
        defC0    = ones(1,numC);
    case 'Sequential'
        Kfun_str = 'diag(-k) + diag(k(1:end-1),-1)'; % Diagonal + diagonal-1
        defC0    = [1 zeros(1, numC-1)];
end

% Convert the K matrix string to function
Kmat_Fun    = str2func(['@(k) ' Kfun_str]);

% Evaluate the K matrix and get the number of species
Ktest       = Kmat_Fun(1:numC);
Ktest       = Ktest(Ktest~=0);
numTaus     = length(unique(abs(Ktest)));

% Check that the K matrix has a consistent size
if numTaus ~= numC
    error('Wrong definition of the K matrix function: K matrix size inconsistent.');
end