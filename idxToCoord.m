function coord=idxToCoord(id,nRows,nCols,im)%x abscisses , y ordonnee
    x = ceil(id/nRows);
    y = mod(id-1,nRows)+1;
    red = double(im(id));
    green = double(im(id+nRows*nCols));
    blue = double(im(id+2*nRows*nCols));
    coord=[x/nCols,y/nRows,red/255,green/255,blue/255];
end