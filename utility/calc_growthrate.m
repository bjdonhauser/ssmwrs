function growthrate_y = calc_growthrate(y, datatype)
% growthrate_y  --  Takes in a set of data and returns the annualized 
%                   growth rate for that data.
%
%****f* SSMWRS/lam_constraints
%
% NAME
%   growthrate_y    --  Takes in a set of data and returns the annualized 
%                       growth rate for that data. 
%
% SYNOPSIS
%   [growthrate_y] = calc_growthrate(y, [datatype])
%
% INPUTS
%   * y             --  data (vector of length N)
%   * datatype      --  string representation of the datatype. Values 
%                       'monthly', 'quarterly', and 'yearly' are 
%                       accepted. (string)      
%
% OUTPUTS
%   * growthrate_y  --  the annualized growth rate of y (vector of length
%                       N-1)
%
% SIDE EFFECTS
%
% DESCRIPTION
%   Takes in a set of data and returns the annualized growth rate for that 
%   data. 
%
% REFERENCES
%
% SEE ALSO
%   lam_main
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-12
%
%***

%=========================================================================%
% Set Defaults
%=========================================================================%

if nargin <= 1 || isempty(datatype)
    datatype = 'yearly';
end

%=========================================================================%
% Calculate
%=========================================================================%

N = length(y);
growthrate_y = log(y(2:N)./y(1:(N-1)));

switch datatype
    
    case {'monthly'}
        growthrate_y = 1200*growthrate_y;
        return
        
    case {'quarterly'}
        growthrate_y = 400*growthrate_y;
        return
        
    case {'yearly'}
        growthrate_y = 100*growthrate_y;
        return
        
end


       
