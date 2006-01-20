function torplotbw(fname)
    set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',0.8,'defaultpatchlinewidth',.7);
    idxDAT = importdata(strcat(fname,'.idx'));
    a = find(idxDAT(1,:) < inf );
    b = find(idxDAT(2,:) < inf );
    idxA = idxDAT(1,a);
    idxB = idxDAT(2,b);
    idxA(length(idxA)+1)=1.0;
    idxB(length(idxB)+1)=1.0;

    sol = importdata(strcat(fname,'.dat'));
    %sol = soltan;
    sol(length(idxB),:) = sol(1,:);
    sol(:,2*length(idxA)-1:2*length(idxA)) = sol(:,1:2);

    Z1 = sol(:,1:2:2*length(idxA));
    Z2 = sol(:,2:2:2*length(idxA));
    X1 = kron(idxA,linspace(1,1,length(idxB))');
    Y1 = kron(idxB,linspace(1,1,length(idxA))')';
    figure(1);
    sf = surf(X1,Z1,Z2,'FaceColor','interp',...
        'EdgeColor','none',...
        'FaceLighting','gouraud');
    shading interp
    colormap(gray);
    axis tight
    ah = xlabel('t/\tau');
    set(ah,'FontSize',16,'FontName','Times');
    ah = ylabel('x(t)/h_0');
    set(ah,'FontSize',16,'FontName','Times');
    ah = zlabel('(dx(t)/dt)/h_0');
    set(ah,'FontSize',16,'FontName','Times');
    figure(2);
    surf(idxA,idxB,Z1);
    return
   