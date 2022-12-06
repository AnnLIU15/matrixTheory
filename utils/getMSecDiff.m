function msec = getMSecDiff(last_time)
%--------------------------------------------------------------------------
%  calculate the sec with msec
%--------------------------------------------------------------------------
%     compute PSNR and reconstruction error for the recovered image and
%     original image
%
%     Inputs:
%         last_time      --- time
%     Outputs:
%         msec           --- sec.msec
%--------------------------------------------------------------------------
    msec = (datenum(datetime('now') ) - datenum(last_time));
    % 24*3600
    msec = msec * 86400;
end