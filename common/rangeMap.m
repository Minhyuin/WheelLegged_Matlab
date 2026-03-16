function [rslt] = rangeMap(scale, min, max)
    period = max - min;
    if (period <= 0)
        rslt = min;
    end
    offset = mod(scale - min, period);
    if (offset < 0)
        offset = offset + period;
    end
    rslt = offset + min;
end