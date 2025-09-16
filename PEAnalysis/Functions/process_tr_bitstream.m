function out = process_tr_bitstream(x,y, options)
%PROCESS_TR_BITSTREAM Summary of this function goes here
%   Detailed explanation goes here

Ts = x(2)-x(1);
y_thr = (max(y)+min(y))/2;  % calculae trigger point threshold
% discriminate over threshold
y_dig = y < y_thr;
bit_l = round(1/options.baud_rate/Ts);
idx=find(y_dig>0.5);
idx_first = idx(1); % find begining of start bit

sampling = (1:options.n_bits)*bit_l+idx_first+round(bit_l/2);
bitstream = y_dig(sampling);

out = binaryVectorToDecimal(bitstream');
end

