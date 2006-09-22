function out = pmat2branch( file )
  load(file);
 
  for i=1:pdde_npoints
    
    %kind can more than 'psol' only.
    branch(i).kind = 'psol';
    branch(i).period = pdde_par(1,i);
    branch(i).parameter = pdde_par(2:end-3,i);
    branch(i).stability.mu = pdde_mul(:,i);
    ndeg = length(pdde_elem(:,i))-1;
    nint = length(pdde_mesh(:,i))-1;
    branch(i).degree = ndeg;
    for k=1:nint
        branch(i).mesh((k-1)*ndeg+1:k*ndeg) = pdde_mesh(k,i).*ones(ndeg,1) + ...
            (pdde_mesh(k+1,i)-pdde_mesh(k,i)).*pdde_elem(1:end-1,i);
    end
    branch(i).mesh(end+1) = pdde_mesh(end,i);
    branch(i).profile = reshape( pdde_prof(:,i), pdde_ndim(i),[]);
  end
  
  %the actual information about the solution is stored in br.point(ii)
  out.point=branch;
  %just as information
  out.method.name='PDDE_CONT';
  
return
