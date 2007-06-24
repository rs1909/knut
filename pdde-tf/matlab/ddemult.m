function mult = ddemult(point)
    mult = [ ];
    for i=1:length(point) 
        mult(i,1:length(point(i).stability.mu))=point(i).stability.mu;
    end

end