function msec = getMSecDiff(last_time)
    msec = (datenum(datetime('now') ) - datenum(last_time));
    % 24*3600
    msec = msec * 86400;
end