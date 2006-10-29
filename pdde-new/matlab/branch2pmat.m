function branch2pmat( file, branch, nmul )
  
  %the actual information about the solution is stored in br.point(ii)
  sol=branch.point;

  pdde_npoints = length( sol );
  pdde_ntrivmul = 1;
  pdde_mul = complex(zeros(nmul, pdde_npoints,'double'));
  for i = 1:length( sol )
    [ndim tmp] = size( sol(i).profile );
    ndeg = sol(i).degree;
    nint = floor( tmp/ndeg );
    pdde_ndim(:,i) = ndim;
    pdde_par(:,i) = [sol(i).period, sol(i).parameter, 0, 0, 0];
    pdde_mul(:,i) = sol(i).stability.mu(1:nmul);
    pdde_elem(:,i) = sol(i).mesh(1:ndeg+1)./sol(i).mesh(ndeg+1);
    pdde_mesh(:,i) = sol(i).mesh(1:ndeg:end);
    pdde_prof(:,i) = reshape( sol(i).profile, 1, ndim*(ndeg*nint + 1) );
  end
  if isreal(pdde_mul)
      pdde_mul = complex(pdde_mul);
  end
  save(file, 'pdde_npoints', 'pdde_ntrivmul', 'pdde_elem', 'pdde_ndim', ...
      'pdde_par', 'pdde_mul', 'pdde_mesh', 'pdde_prof', '-V4')
return
