function [Kmat_Fun,numTaus,numC,defTbl,defRowNames] = ParseKmatrix(modelIdx)

% This function defines a kinetic target model for first-order coupled reaction systems.
% A few models have been pre-defined, and a custom model can be typed.
%
% To correctly display the reaction graph, all rates must be written following the order of the species:
% e.g. Write all rates that involve species A, then species B, and so on. This can be counter-intuitive sometimes.
%
% Ricardo Fernández-Terán /  2022.05.01 / v1.3b


% Define the K matrix from target model
switch modelIdx
    case 1  % 'A <=> B';
        Kfun_str = '[-k(1)   k(2); k(1)  -k(2)]';
    case 2  % 'A <=> B; A -> ; B -> ';
        Kfun_str = '[-(k(1)+k(2)) k(3); k(1) -(k(3)+k(4))]';
    case 3  % 'A <=> B; A ->';         
        Kfun_str = '[-(k(1)+k(2)) k(3); k(1) -k(3)]';
    case 4  % 'A <=> B; B ->';         
        Kfun_str = '[-k(1) k(2); k(1) -(k(2)+k(3))]';
    case 5  % 'A <=> B; B -> C; C -> ';   
        Kfun_str = '[-k(1) k(2) 0 ; k(1) -(k(2)+k(3)) 0; 0 k(3) -k(4)]';
    case 6  % 'A <=> B; B<=>C; C ->';
        Kfun_str = '[-k(1) k(2) 0; k(1) -(k(2)+k(3)) k(4); 0 k(3) -(k(4)+k(5))]'; 
    case 7  % 'A -> B; B<=>C; B -> ; C ->';
        Kfun_str = '[-k(1) 0  0; k(1) -(k(2)+k(3)) k(4); 0 k(2) -(k(4)+k(5))]';
    case 8  % 'A -> B -> D; A -> C -> D; D ->';
        Kfun_str = '[-k(1)-k(2) 0 0 0; k(2) -k(4) 0 0; k(1) 0 -k(3) 0; 0 k(4) k(3) -k(5)]';
    case 9  % 'A -> B -> C ; A -> D -> E';
        Kfun_str = '[-k(1)-k(2) 0 0 0 0; k(1) -k(3) 0 0 0; 0 k(3) 0 0 0; k(2) 0 0 -k(4) 0; 0 0 0 k(4) 0]';
    case 10 % 'A -> B -> C -> ; A -> D -> E ->';
        Kfun_str = '[-k(1)-k(2) 0 0 0 0; k(1) -k(3) 0 0 0; 0 k(3) -k(4) 0 0; k(2) 0 0 -k(5) 0; 0 0 0 k(5) -k(6)]';
    case 11 % 'A -> B -> ; A -> C ->'
        Kfun_str = '[-k(1)-k(2) 0 0; k(1) -k(3) 0; k(2) 0 -k(4)]';
    case 12 % 'A -> B -> C <=> D; C ->; D ->'
        Kfun_str = '[-k(1), 0, 0, 0; k(1), -k(2), 0, 0; 0, k(2), -(k(3) + k(4)), k(5); 0, 0, k(3), -(k(5) + k(6))]';
    case 13 % Custom model
        prompt   = 'Enter K matrix in functional form [i.e. as a function of a vector "k" containing the rate constants (see example)]:';
        title    = 'Define K matrix for target model';
        dims     = [10 50];
        definput = {['[-k(1)-k(2) k(3)'; 'k(1) -k(2)-k(3)]']};

        Kfun_inp = inputdlg(prompt,title,dims,definput);
        
        % Parse the input (a char array) into a properly formatted function string
        Kfun_str = Kfun_inp{1}(1,:);
        for i=2:size(Kfun_inp{1},1)
            Kfun_str = [Kfun_str ';' Kfun_inp{1}(i,:)]; %#ok<*AGROW> 
        end
end

% Check how many parameters (rate constants) are needed
matchStr    = regexp(Kfun_str,'k\(\d\)','match');
numTaus     = numel(unique(matchStr));

% Convert the K matrix string to function
Kmat_Fun    = str2func(['@(k) ' Kfun_str]);

% Evaluate the K matrix and get the number of species
numC        = unique(size(Kmat_Fun(ones(numTaus,1))));

% Check that the K matrix is square
if length(numC) > 1
    error('Wrong definition of the K matrix function: the K matrix is not square!');
end

% Generate new default rate table (for the GUI)
defTbl      = [];
defRowNames = {};
defKin      = [1 Inf 0];

for i=1:numTaus
    defTbl      = [defTbl; [num2cell(defKin*5*i) {false}]];
    defRowNames = [defRowNames; {num2str(i)}]; 
end
