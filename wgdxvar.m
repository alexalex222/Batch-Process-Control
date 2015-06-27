function s = wgdxvar(varargin)

% It is a function to transfer/convert variables

% The first input is the name of the variable
s.name = varargin{1};

% The second input is the value of the variable
s.val = varargin{2};

% Variable type
s.type = 'parameter';

% Variable form
s.form = 'full';

% The third input is dimension of the variable
s.dim = varargin{3};

end