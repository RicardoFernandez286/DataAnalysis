function VisualiseKMat(Kmat_Fun,numTaus,numC,modelType)

% [Kmat_Fun,numTaus,numC] = ParseKmatrix(1);
AMat = Kmat_Fun(1:numTaus);

% LMat = Kmat_Fun(1:numTaus);
% syms p(numTaus)
% sum(Kmat_Fun(p),1)
% AMat = diag([1,1,1,1],-1)-eye(5);  numC = 5; numTaus=5;

% if not mass conservation, then add an additional row and column with the difference
nodenames = string;
for i=1:numC
    nodenames = [nodenames char(64+i)];
end

if sum(AMat,'all') ~= 0
    A2  = AMat - diag(diag(AMat));
    BMat= [[A2; -sum(AMat,1)] zeros(numC+1,1)];
    nodenames = [nodenames '0'];
else
    BMat= AMat - diag(diag(AMat));
end
nodenames = nodenames(2:end);

% edgelabels = {'k_{1}','k_{3}','k_{2}','k_{4}'};
% edgelabels = {'k_{1}','k_{3}','k_{2}','k_{4}'};

edgelabels = string;
for i=1:numTaus
    edgelabels = [edgelabels string(['k_{' num2str(i) '}'])];
end
edgelabels = edgelabels(2:end);

G = digraph(BMat',nodenames,'omitselfloops');

fh=figure(1);
clf(fh);
fh.Name = 'Kinetic Network Graph';
movegui(fh,'west')
fh.Color='w';
fh.Position(3:4) = [390 420];

w = [1 1 1];

switch modelType
    case {'Parallel','Sequential'}
        if numTaus == 2
            A = [0 0 0.85; 0.85 0 0];
        else
            A = brighten(turbo(numTaus+1),-0.4);
            A = A(1:end-1,:);
        end
        
        edgeLabelCol = A;
        edgeCol = A;
        nodeCol = 0.75*w;
    case 'Target'
        if numC == 2
            A = [0 0 0.85; 0.85 0 0];
        else
            A = brighten(turbo(numC+1),-0.4);
            A = A(1:end-1,:);
        end
        
        if sum(AMat,'all') ~= 0
                A = [A; 0.85 0.85 0.85];
        end
        
        edgeLabelCol = 0*w;
        edgeCol = 0.75*w;
        nodeCol = A;
    otherwise
        return
end

plot(G,'-o','NodeColor',nodeCol,'Layout','layered',...
    'EdgeLabel',edgelabels,'LineWidth',2, ...
    'NodeFontSize',14,'NodeFontWeight','bold', 'EdgeColor',edgeCol,...
    'EdgeFontSize',14,'EdgeFontWeight','bold','EdgeLabelColor',edgeLabelCol, ...
    'ArrowSize',15,'ArrowPosition',0.75,'MarkerSize',10,'EdgeAlpha',1)

