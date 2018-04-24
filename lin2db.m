function varargout = lin2db(varargin)

for i = 1:nargin
    varargout{i} = 10*log10(varargin{i});
end

end

