fprintf('Path scores of subway networks by city. \n')

file1='Shanghai-2010-adjacency.txt';
C=adjacency_matrix(file1);   
disp('path scores for the stations in 2010 Shanghai')
nonZero(C);
graphProperties(C);
export(file1, 'shanghaiPathName.xlsx')

file2 = 'Osaka-2010-adjacency.txt';
D = adjacency_matrix(file2);
disp('path scores for the stations in 2010 Osaka')
nonZero(D);
graphProperties(D);
export(file2, 'osakaPathName.xlsx')

file3='Tokyo-2009-adjacency.txt';
E = adjacency_matrix(file3);
disp('path scores for the stations in 2009 Tokyo')
nonZero(E)
graphProperties(E);
export(file3, 'tokyoPathName.xlsx')

file4 = 'Seoul-2009-adjacency.txt';
F = adjacency_matrix(file4);
disp('path scores for the stations in 2009 Seoul')
nonZero(F);
graphProperties(F);
export(file4, 'seoulPathName.xlsx')

file5 = 'Moscow-2009-adjacency.txt';
G = adjacency_matrix(file5);
disp('path scores for the stations in 2009 Moscow')
nonZero(G);
graphProperties(G);
export(file5, 'moscowPathName.xlsx')

file6 = 'Barcelona-2010-adjacency.txt';
H = adjacency_matrix(file6);
disp('path scores for the stations in 2010 Barcelona')
nonZero(H);
graphProperties(H);
export(file6, 'barcelonaPathName.xlsx')

file7 = 'NYC-2009-adjacency.txt';
I = adjacency_matrix(file7);
disp('path scores for the stations in 2009 NYC');
nonZero(I);
graphProperties(I);
export(file7, 'nycPathName.xlsx')

file8 = 'Berlin-2010-adjacency.txt';
J = adjacency_matrix(file8);
disp('path scores for the stations in 2010 Berlin');
nonZero(J);
graphProperties(J);
export(file8, 'berlinPathName.xlsx')

file9 = 'HongKong-2009-adjacency.txt';
K = adjacency_matrix(file9);
disp('path scores for the stations in 2009 Hong Kong');
nonZero(K);
graphProperties(K);
export(file9, 'hongkongPathName.xlsx')

file10 = 'Mexico-2009-adjacency.txt';
L = adjacency_matrix(file10);
disp('path scores for the stations in 2009 Mexico');
nonZero(L);
graphProperties(L);
export(file10, 'mexicoPathName.xlsx')

file11 = 'London-2009-adjacency.txt';
M = adjacency_matrix(file11);
disp('path scores for the stations in 2009 London');
nonZero(M);
graphProperties(M);
export(file11, 'londonPathName.xlsx')

file12 = 'Paris-2009-adjacency.txt';
N = adjacency_matrix(file12);
disp('path scores for the stations in 2009 Paris');
nonZero(N);
graphProperties(N);
export(file12, 'parisPathName.xlsx')

file13 = 'Beijing-2010-adjacency.txt';
O = adjacency_matrix(file13);
disp('path scores for the stations in 2010 Beijing');
nonZero(O);
graphProperties(O);
export(file13, 'beijingPathName.xlsx')

file14 = 'Madrid-2009-adjacency.txt';
P = adjacency_matrix(file14);
disp('path scores for the stations in 2009 Madrid');
nonZero(P);
graphProperties(P);
export(file14, 'madridPathName.xlsx')

function export(dataName, filename)
    A = adjacency_matrix(dataName);

    fileID = fopen(dataName);
    S = textscan(fileID,'%s %s');
    fclose(fileID);

    G = digraph(transpose(S{2}),transpose(S{1}));

    pathScore = pathscore_SHL(A);
    n=numnodes(G);
    node = string.empty(n,0);
    
    for i=1:n
        node(i)=G.Nodes.Name(i);
    end

    T = table(transpose(pathScore),transpose(node));
    writetable(T, filename, 'Sheet','Range');
end

% input: adjacency matrix
% output: list of non-zero path scores and their index(station number)
function nonZero(A)
    B = pathscore_SHL(A);    % find path scores of all nodes
    
    fprintf('\n')

    % print out the nodes that have non-zero path scores
    total = 0;
    for i = 1:length(B)
        if B(i) > 0
            %fprintf('node %d : %f\n',i,B(i));
            total = total + 1;
        end
    end

    fprintf('There are %d nodes with nonzero path score \n', total)
    fprintf('The highest path score is %f', max(B))
end


% input: .txt files for subway adjacency numbers, AND number of subway stations
% output: adjacency matrix of 0's and 1's 
function [A] = adjacency_matrix(G)
    filename= G;
    fileID = fopen(filename);

    % separates the text file into a source vector S{1} and destination
    % vector S{2}
    S = textscan(fileID,'%s %s');
    fclose(fileID);
    
    G = digraph(transpose(S{2}),transpose(S{1}));

    A = adjacency(G);
end

function [nodes, edges, avgdeg, bighub, diam, avgpath, dense] = graphProperties(G);
%input is adjacency matrix G, output is a lot of network properties
H = graph(G);
nodes = numnodes(H);
edges = numedges(H);
avgdeg = (2*edges)/nodes;
bighub = max(degree(H));
d = distances(H);
diam = max(d(:));
avgpath = mean(d(:));
dense = nnz(adjacency(H))./numel(adjacency(H));
fprintf('\nGraph Properties\n');
fprintf('\tNumber of Nodes: %f\n',nodes);
fprintf('\tNumber of Edges: %f\n',edges);
fprintf('\tMean Degree: %6.4f\n',avgdeg);
fprintf('\tDegree of Largest Station: %6.4f\n',bighub);
fprintf('\tDiameter: %6.4f\n',diam);
fprintf('\tMean Path Length: %6.4f\n',avgpath);
fprintf('\tDensity: %6.4f\n',dense);
fprintf('\n');
end


function [Y] = pathscore_SHL (G)
n=length(G);
Y(1:n)=0;

    function [d,Q]=dist1(v,w,d,Q)
        if d(w)<0
            Q=[Q w];
            d(w)=d(v)+1;
        end
    end
    function sig=pathcount(v,w,d,sig)
        if d(w)<d(v)
            sig(w)=sig(w)+sig(v);
        end
    end

%for s=1:n
%    for t=1:(s-1)
%        if G(s,t)>0
[ex,ey] = find(triu(G) > 0);
EDGES = [ex ey];
m = length(EDGES);
for e=1:m
    	if mod(e,round(m/10)) == 0    
            fprintf('');    end   
    	    s=EDGES(e,1); t=EDGES(e,2);
            a=G(s,t);
            G(s,t)=0;
            G(t,s)=0;
            sigs(1:n)=0;
            sigt(1:n)=0;
            sigs(s)=1;
            sigt(t)=1;
            ds(1:n)=-1;
            dt(1:n)=-1;
            ds(s)=0;
            dt(t)=0;
            Q=s;
            while ~isempty(Q)
               v=Q(1);
               Q(1)=[];
               for w=find(G(v,:))
                   [ds,Q]=dist1(v,w,ds,Q);
               end
            end
            Q=t;
            while ~isempty(Q)
               v=Q(1);
               Q(1)=[];
               for w=find(G(v,:))
                   [dt,Q]=dist1(v,w,dt,Q);
                   sigt=pathcount(v,w,ds,sigt);
               end
            end
            Q=s;
            while ~isempty(Q)
               v=Q(1);
               Q(1)=[];
               for w=find(G(v,:))
                   if dt(w)<dt(v)
                       Q=[Q w];
                   end
                   sigs=pathcount(v,w,dt,sigs);
               end
            end
            for v=setdiff(1:n,[s t])
               if(sigs(v)*sigt(v) > 0)
                 Y(v)=Y(v)+sigs(v)*sigt(v)/sigs(t);
               end
            end
            G(s,t)=a;
            G(t,s)=a;
%        end
%    end
end

end
