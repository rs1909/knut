idxDAT = importdata('soltan.idx');
a = find(idxDAT(1,:) < inf );
b = find(idxDAT(2,:) < inf );
idxA = idxDAT(1,a);
idxB = idxDAT(2,b);
load soltan.dat
sol = soltan;
Z1 = sol(:,1:2:2*length(idxA));
Z2 = sol(:,2:2:2*length(idxA));
X1 = kron(idxA,linspace(1,1,length(idxB))');
Y1 = kron(idxB,linspace(1,1,length(idxA))')';
figure(1);
surf(X1,Z1,Z2,'FaceColor','interp',... %    'EdgeColor','none',...
    'FaceLighting','phong');
camlight left
figure(2);
surf(idxA,idxB,Z1,'FaceColor','interp',... %    'EdgeColor','none',...
    'FaceLighting','phong');
camlight left