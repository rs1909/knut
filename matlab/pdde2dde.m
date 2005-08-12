function branch = pdde2dde( file )
    in = load( file );
    [nsol, dummy] = size( in );
    for i=1:nsol
        npar = in(i,1);
        branch(i).kind = 'psol';
        branch(i).period = in(i,2);
        branch(i).parameter = in(i,3:1+npar);
        nmul = in(i,npar+2);
        branch(i).stability.mu = (in(i,npar+3:2:npar+2*nmul+2) ...
                             + j*in(i,npar+4:2:npar+2*nmul+2))';
        ndim = in(i,npar+2*nmul+3);
        nint = in(i,npar+2*nmul+4);
        ndeg = in(i,npar+2*nmul+5);
        branch(i).mesh = in(i,npar+2*nmul+6:npar+2*nmul+5+ndeg*nint+1);
        branch(i).profile = reshape( in(i,npar+2*nmul+6+ndeg*nint+1:end), ...
           ndim, ndeg*nint+1);
    end
end