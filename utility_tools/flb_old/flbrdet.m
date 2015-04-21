function [num_chan, chan_len, chan_format, domaine_name, ...
	incr, start, comment, ...
	chan_name, chan_min, chan_max ] = flbrdet (fid)
% read details from a flb file
[a,count]=fread(fid,4,'int8');
if (count ~= 4),
	num_chan=-1;
	disp('EOF');
	return
end	
flb_version = iaxp2mac(a);
% if (flb_version ~= 2),
% 	disp ('Header error')
% 	num_chan=-2;
% 	return
% end
num_chan = iaxp2mac(fread(fid,4,'int8'))
num_real = iaxp2mac(fread(fid,4,'int8'))
chan_len = iaxp2mac(fread(fid,4,'int8'))
chan_format = iaxp2mac(fread(fid,4,'int8'))


dlen = iaxp2mac(fread(fid,4,'int8'))
domaine_name = setstr(fread (fid, dlen,'uchar'))
incr = fread(fid,1,'float32')
start = fread(fid,1','float32')
clen = iaxp2mac(fread(fid,4,'int8'))
comment = setstr(fread(fid,clen,'uchar'))
chan_name=[];
for i=1:num_chan
    [inlen, count]=fread(fid,4,'int8')
    nlen = iaxp2mac(inlen)
  chan_name = [chan_name setstr(fread(fid,nlen,'char'))' '|']
end
chan_min =fread (fid,num_chan(1),'float64')
chan_max =fread(fid,num_chan(1),'float64')
return







