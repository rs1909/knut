function torplot(fname,n)
    load(fname)
    if knut_npoints+1 < n
        disp 'Error: The specified point does not exist in the file'
        return
    end
    idxA = knut_mesh1(:,n);
    idxB = knut_mesh2(:,n);
    idxA(end+1)=1.0;
    idxB(end+1)=1.0;
    
    ndim  = knut_ndim(n);
    nint1 = knut_nint1(n);
    nint2 = knut_nint2(n);
    ndeg1 = knut_ndeg1(n);
    ndeg2 = knut_ndeg2(n);
%    when already reshaped in the software
%    sol_t = knut_blanket(:,n);
%    when it is not reshaped in the software
    sol_t = reshape(knut_blanket(:,n),ndim,ndeg1,ndeg2,nint1,nint2);
    sol = reshape(permute(sol_t,[2,4,3,5,1]),ndeg1*nint1,ndeg2*nint2,ndim);
    sol(end+1,:,:) = sol(1,:,:);
    sol(:,end+1,:) = sol(:,1,:);
    figure(1);
    surf(idxA,idxB,squeeze(sol(:,:,1)));
    figure(2);
    surf(idxA,idxB,squeeze(sol(:,:,2)));
    
    Z1 = squeeze(sol(:,:,1));
    Z2 = squeeze(sol(:,:,2));
    X1 = kron(idxA,linspace(1,1,length(idxB)));
    Y1 = kron(idxB,linspace(1,1,length(idxA))')';
    figure(3);
    surf(X1,Z1,Z2,'FaceColor','interp',...    
        'EdgeColor','none',...
        'FaceLighting','gouraud');
end
