function p = paramlog2array(s)
%Transforms general parameter log returned from *load_logfile.m* to a
%numerical array.
%
%   s: struct-array with fields
%       id
%       datetime
%       params
%
%   p: double-array
%
%
%   NOTE: non-numeric params are mapped to an NaN.

nshot=numel(s);                 % number of shots logged
nparam=numel(s(1).params);      % length of parameter vector

p=NaN(nshot,1+nparam);      % initialise
p(:,1)=[s.id]';             % shot ID
% convert string param vector to num-vect
for ii=1:nshot
    for jj=1:nparam
        p(ii,1+jj)=str2double(s(ii).params{jj});
    end
end

end