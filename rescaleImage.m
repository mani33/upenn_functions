function imageMatrix = rescaleImage( imageMatrix, varargin )
%imageMatrix = rescaleImage( imageMatrix,'qval', quantVals )
%   rescaleImage( imageMatrix,'qval', [0.01 0.99] )
args.qval =  [0.01 0.99];
args = parseVarArgs(args,varargin{:});

co = class(imageMatrix);
if isa(co,'double')
    ma = 1;
else
    ma = double(intmax(co));
end

im = double(imageMatrix);
q = quantile(im(:), args.qval);
im = im-q(1);
im = (im/q(2)) * ma;%#ok
imageMatrix = eval(sprintf('%s(im)',co));


