function branch2pmat( file, branch, nmul )
  
  %the actual information about the solution is stored in br.point(ii)
  sol=branch.point;

  knut_npoints = length( sol );
  knut_ntrivmul = 1;
  knut_mul = complex(zeros(nmul, knut_npoints,'double'));
  for i = 1:length( sol )
    [ndim tmp] = size( sol(i).profile );
    ndeg = sol(i).degree;
    nint = floor( tmp/ndeg );
    knut_ndim(:,i) = ndim;
    knut_par(:,i) = [sol(i).period, sol(i).parameter, 0, 0, 0];
    knut_mul(:,i) = sol(i).stability.mu(1:nmul);
    knut_elem(:,i) = sol(i).mesh(1:ndeg+1)./sol(i).mesh(ndeg+1);
    knut_mesh(:,i) = sol(i).mesh(1:ndeg:end);
    knut_prof(:,i) = reshape( sol(i).profile, 1, ndim*(ndeg*nint + 1) );
  end
  if isreal(knut_mul)
      knut_mul = complex(knut_mul);
  end
  save(file, 'knut_npoints', 'knut_ntrivmul', 'knut_elem', 'knut_ndim', ...
      'knut_par', 'knut_mul', 'knut_mesh', 'knut_prof', '-V4')
return
