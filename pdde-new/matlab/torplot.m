function torplot(fname,n)
    load(fname)
    if pdde_npoints+1 < n
        disp 'Error: The specified point does not exist in the file'
        return
    end
    idxA = pdde_mesh1(:,n);
    idxB = pdde_mesh2(:,n);
    idxA(end+1)=1.0;
    idxB(end+1)=1.0;
    
    ndim  = pdde_ndim(n);
    nint1 = pdde_nint1(n);
    nint2 = pdde_nint2(n);
    ndeg1 = pdde_ndeg1(n);
    ndeg2 = pdde_ndeg2(n);
%    when already reshaped in the software
%    sol_t = pdde_blanket(:,n);
%    when it is not reshaped in the software
    sol_t = reshape(pdde_blanket(:,n),ndim,ndeg1,ndeg2,nint1,nint2);
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
