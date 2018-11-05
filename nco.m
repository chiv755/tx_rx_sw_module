function out_nco = nco(frequency)
%Initialize the sine table    
	i = 0:255;
    sintable32 = round(127*(sin(2*pi*i/256))+128*ones(size(i))+0.5);
    
    accumulator = 0;
    value = sintable32(1);
   
    %256 table entries, 24 bits of fractional index; 2^24 = 16777216
	phase = round(frequency*256*16777216+0.5);
    
    for i = 0:200 
			index = (accumulator >> 24) & 0xff;
            value = sintable32(index);
            accumulator = accumulator + phase;
    end
end