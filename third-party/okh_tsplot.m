function okh_tsplot(y, info, varargin)
%---------------------------------------------------------------------%
% PURPOSE: time-series plot with dates
%---------------------------------------------------------------------%
% USAGE:  okh_tsplot(y, info, varargin) 
%           where: y      = vector of series to be plotted
%                  info   = plot information
%                   * info.yy   = the beginning year
%                   * info.pr   = the beginning period of the year
%                   * info.freq = frequency : Y, Q, M
%                     info.dateform  = dateform (see datetick)
%                     info.fsize  = font size
%                     info.hline  = horizontal line
%                     info.hlineColor  = horizontal line color
%                     info.bcycle = Draw business cycle (1, [0])
%                     info.bcycleColor = Color specification for business cycle recession 
%                  varargin = line characteristics
%            '*' means requirements.
% AUTHOR:
%  Kum Hwa Oh
%---------------------------------------------------------------------%


%---------------------------------------------------------------------%
% Initialization
%---------------------------------------------------------------------%
[obs1, obs2] = size(y) ;
if obs1 > obs2
    n_obs = obs1 ;
else
    n_obs = obs2 ;
end

if ~isfield(info, 'fsize')
    info.fsize = 9 ;
end
if ~isfield(info, 'dateform')
    if info.freq == 'Y' 
        info.dateform = 'yyyy' ;
    elseif info.freq == 'Q'
        info.dateform = 'QQ-YY' ;
    else
        info.dateform = 'mmmyy'
    end
end

if ~isfield(info, 'bcycle')
    info.bcycle = 0 ;
end
if ~isfield(info, 'bcycleColor')
    info.bcycleColor = [0.4941 0.4941 0.4941] ;
end

recession = [1948 11 1949 10 ;
             1953 07 1954 05 ;
             1957 08 1958 04 ;
             1960 04 1961 02 ;
             1969 12 1970 11 ; 
             1973 11 1975 03 ;
             1980 01 1980 07 ;
             1981 07 1982 11 ;
             1990 07 1991 03 ;
             2001 03 2001 11 ] ;



%---------------------------------------------------------------------%
% make dates
%---------------------------------------------------------------------%
switch info.freq 
    
    case {'Y'}  % case of annual series 
        yr  = info.yy ;
        for i = 1:n_obs
            yrs(i,1) = yr ;
            mths(i,1) = 1 ;
            yr = yr + 1 ;
        end;
   
     case {'Q'}  % case of quarterly series
        yr  = info.yy ;
        if info.pr == 1 
            mth = 1 ;
        elseif info.pr == 2
            mth = 4 ;
        elseif info.pr == 3
            mth = 7 ;
        else
            mth = 10 ;
        end;

        for i = 1:n_obs
            yrs(i,1) = yr ;
            mths(i,1) = mth ;
            mth = mth + 3 ;
            if mth > 12
                yr = yr + 1 ;
                mth = 1;
            end;
        end;
            

     case {'M'}  % case of monthly series
        yr  = info.yy ;
        mth = info.pr ;
        for i = 1:n_obs
            yrs(i,1) = yr ;
            mths(i,1) = mth ;
            mth = mth + 1 ;
            if mth > 12
                yr = yr + 1 ;
                mth = 1;
            end;
        end;

end


%---------------------------------------------------------------------%
% plot the series
%---------------------------------------------------------------------%
plot(datenum(yrs,mths,1), y, varargin{:}) ;
hold on

% draw horizontal line
if isfield(info, 'hline') & isnumeric(info.hline)
    if ~isfield(info, 'hlineColor')
        info.hlineColor = 'r' ;
    end
    plot(datenum(yrs,mths,1), ones(1,n_obs)*info.hline, info.hlineColor) 
end


%---------------------------------------------------------------------%
% shade recession periods
%---------------------------------------------------------------------%
if info.bcycle == 1
    limY = axis ;
    if ~isfield(info, 'bcycleY') 
        info.bcycleY = [limY(1,3), limY(1,4)] ;
    end
    if ischar(info.bcycleY)
        info.bcycleY = [limY(1,3), limY(1,4)] ;
    end

    i = 1 ;
    j = 1 ;
    flag_s = 0 ;
    flag_e = 0 ;
    
    for j = 1:size(recession,1)
        recess_st = datenum(recession(j,1),recession(j,2),1) ;
        recess_ed = datenum(recession(j,3),recession(j,4),1) ;
        if info.freq == 'Q'
            recess_st = recess_st - 80 ;    
        end
        
        if (datenum(yrs(i,1),mths(i,1),1) <= recess_ed)
            k = i ;
            for i = k:n_obs-1
                if (datenum(yrs(i,1),mths(i,1),1) >= recess_st)*(flag_s == 0)
                    recess_start = datenum(yrs(i,1),mths(i,1),1) ;
                    flag_s = 1 ;
                end
                if (datenum(yrs(i+1,1),mths(i+1,1),1) > recess_ed)
                    recess_end = datenum(yrs(i,1),mths(i,1),1) ;
                    flag_e = 1 ;
                    break
                end
            end    
        end

        if flag_e == 1
            %datestr([recess_start, recess_end], 27)
            fill([recess_start, recess_start, recess_end, recess_end], ...
                 [info.bcycleY(1,2), info.bcycleY(1,1), info.bcycleY(1,1), info.bcycleY(1,2)], ...
                 info.bcycleColor, ...
                 'EdgeAlpha', 0.5, ...
                 'FaceAlpha', 0.5, ...
                 'EdgeColor', info.bcycleColor) ;
            flag_s = 0 ;
            flag_e = 0 ;
        end
    end
end
        
    

%---------------------------------------------------------------------%
% graph settings
%---------------------------------------------------------------------%

set(gca,'fontsize', info.fsize)
set(gca,'tickdir', 'in') 
datetick('x', info.dateform) 
set(gca, 'xcolor', 'blue') 
set(gcf, 'Color', 'w')

hold off


return