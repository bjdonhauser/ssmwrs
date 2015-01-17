function [results_ts] = extract_results_ts(results_element, subelement, start)
% extract_results_ts    --  A utility function that extracts the time 
%                           series of results_element (which is a field of 
%                           the 'results' structure returned by 
%                           ms_ssm_optimizer) as denoted by subelement. 
%
%
%****f* SSMWRS/utility/extract_results_ts
%
% NAME
%   extract_results_ts  --  A utility function that extracts the time 
%                           series of results_element (which is a field 
%                           of the 'results' structure returned by 
%                           ms_ssm_optimizer) as denoted by subelement. 
%
% SYNOPSIS
%   [results_ts] = extract_results_ts(results_element, [subelement], [start])
%
% INPUTS
%   * results_element   --  A field of the a "results" structure returned
%                           by ms_ssm_optimizer
%   * subelement        --  (optional) When the field refers to "S", 
%                           subelement refers to the state. When the field 
%                           refers to anything else, subelement refers to 
%                           the element of a representative vector.
%                           Default: subelement=1.
%   * start             --  (optional). Refers to the starting point of the
%                           time series vector to be returned. This may be
%                           desirable if the first 1:start elements
%                           returned by the kim_filter are unreliable
%                           Default: start=1;
%
% OUTPUTS
%   * results_ts        --  the resulting time series (vector of length
%                           N-start+1)
%
% SIDE EFFECTS
%
% DESCRIPTION
%   A utility function that extracts the time series of results_element 
%   (which is a field of the 'results' structure returned by 
%   ms_ssm_optimizer) as denoted by subelement.
%   
% REFERENCES
%   
%
% SEE ALSO
%   ms_ssm_optimizer, kim_filter
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-12
%
%***

%=========================================================================%
% Set defaults
%=========================================================================%
if nargin <= 2 || isempty(start)
    start=1;
end

if nargin <= 1 || isempty(subelement)
    subelement=1;
end

%=========================================================================%
% Initialize output
%=========================================================================%
N = length(results_element);
results_ts = zeros((N-start+1),1);

%=========================================================================%
% Calculations
%=========================================================================%
% Case in which we're dealing with an "S" field of results. In this case we
% have to use different methods to extract the time series.
if iscell(results_element{1})
    for i=start:N
        results_ts((i-start+1),1) = results_element{i}{subelement};
    end

    % Case in which we're dealing with any field of results besides "S"
else
    for i=start:N
        results_ts((i-start+1),1) = results_element{i}(subelement);
    end
end
