function out = pmat2branch( file )
  load(file);
 
  for i=1:knut_npoints
    
    %kind can more than 'psol' only.
    branch(i).kind = 'psol';
    branch(i).period = knut_par(1,i);
    branch(i).parameter = knut_par(2:end-3,i);
    branch(i).stability.mu = knut_mul(:,i);
    ndeg = length(knut_elem(:,i))-1;
    nint = length(knut_mesh(:,i))-1;
    branch(i).degree = ndeg;
    for k=1:nint
        branch(i).mesh((k-1)*ndeg+1:k*ndeg) = knut_mesh(k,i).*ones(ndeg,1) + ...
            (knut_mesh(k+1,i)-knut_mesh(k,i)).*knut_elem(1:end-1,i);
    end
    branch(i).mesh(end+1) = knut_mesh(end,i);
    branch(i).profile = reshape( knut_prof(:,i), knut_ndim(i),[]);
  end
  
  %the actual information about the solution is stored in br.point(ii)
  out.point=branch;
  %just as information
  out.method.name='KNUT_CONT';
  
return
