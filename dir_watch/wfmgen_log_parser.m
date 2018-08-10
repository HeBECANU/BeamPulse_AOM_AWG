function [uniq_params,id_par,shotid,Ipar] = wfmgen_log_parser(path_to_logfile)
% Reads the log-file created by scanning with wfmgen script
% DKS
% 20180810
%

% load csv as array
%   col1 is shotID; col2.. are params-vector
param_log=load_logfile(path_to_logfile);
param_array=paramlog2array(param_log);


% get unique param-vecs and tag each shot with param-ID
[uniq_params,~,Ipar] = unique(param_array(:,2:end),'rows');
shotid=param_array(:,1);
nparam=size(uniq_params,1);      % number of unique param-set


% group shot-ids by exp-param
id_par=cell(1,nparam);
for ii=1:nparam
    id_par{ii}=shotid(Ipar==ii);
end

end