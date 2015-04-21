function x = strcvrt ( s, form );
if length(s) ==0,
        x=nan;
    elseif isnan(s),
        x=nan;
else
        x=sscanf(s,form);
    end
