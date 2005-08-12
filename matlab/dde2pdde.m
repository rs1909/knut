function dde2pdde( file, sol, nmul )
    fp = fopen( file, 'W' );
    for i = 1:length( sol )-1
        fprintf( fp, '%i\t', length( sol(i).parameter ) + 1 );
        fprintf( fp, '%.16e\t', sol(i).period );
        for j = 1:length( sol(i).parameter )
            fprintf( fp, '%.12e\t', sol(i).parameter(j) );
        end
        fprintf( fp, '%i\t', nmul );
        mmu = [];
        mmu(1:2:2*nmul) = ...
            [ real( sol(i).stability.mu ); ...
            zeros(nmul-length(sol(i).stability.mu),1) ];
        mmu(2:2:2*nmul) = ...
            [ imag( sol(i).stability.mu ); ...
            zeros(nmul-length(sol(i).stability.mu),1) ];
        fprintf( fp, '%.16e\t', mmu );
        [ndim tmp] = size( sol(i).profile );
        ndeg = sol(i).degree;
        nint = floor( tmp/ndeg );
        if ndeg*nint + 1 ~= tmp 
            print 'error'
        end
        fprintf( fp, '%i\t', ndim );
        fprintf( fp, '%i\t', nint );
        fprintf( fp, '%i\t', ndeg );
        fprintf( fp, '%.16e\t', sol(i).mesh );
        flatsol = reshape( sol(i).profile, 1, ndim*(ndeg*nint + 1) );
        fprintf( fp, '%.16e\t', flatsol );
        fprintf( fp, '\n' );
    end
    fclose( fp );
end