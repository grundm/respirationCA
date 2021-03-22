function resp_freq = count_resp_cond(cond_seq,resp)
% resp_freq = count_resp(cond_seq,resp) counts the number responses codes
% (resp) per condition (cond_seq) and compute their frequency
%
% Input:
% cond_seq  - vector with condition sequence (e.g., nt_data.seq(:,2))
% resp      - vector of response sequence (e.g., nt_data.resp1)
%
% Ouput matrix:
% col1      - condition
% col2      - mean of response codes
% col3-x    - frequency of each response code
% colx-y    - number of each response code
% col_end   - total number of responses
%
% Author:           Martin Grund
% Last update:      December 1, 2015

conditions = unique(cond_seq);
resp_freq = zeros(length(conditions),3+2*numel(unique(resp)));

for i = 1:length(conditions)
    resp_tmp = resp(cond_seq==conditions(i));
    tbl_tmp = tabulate(resp_tmp); % not all conditions have all response codes
    resp_freq(i,1:end) = [conditions(i) mean(resp_tmp) tbl_tmp(:,3)'/100 tbl_tmp(:,2)' numel(resp_tmp)]
end