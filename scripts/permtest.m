function [h, p, diff, diff_null] = permtest(a,varargin)
%   Non-parametric permutation test to compare means using
%   within-participant shuffling. Use as a paired-sample permutation test. 
%
%   H0: the two populations have equal means
%   HA: the two populations have unequal means
%
%   Usage:
%     [h p diff diff_null] = permtest(a)                          (compare to zero)
%     [h p diff diff_null] = permtest(a,b)                        (compare to a particular value, or compare the mean of two distributions)
%     [h p diff diff_null] = permtest(a,b,npermutes,pthresh,tail,dim)
%
%   Input:
%     a:          distribution 1
%     b:          (optional) distribution 2, or a single value. default: 0
%     npermutes:  (optional) the number of permutations. default: 1000
%     pthresh:    (optional) the p-threshold for significance. default: 0.05
%     tail:       (optional) test one-tailed (specify 'right' or 'left') or two-tailed (specify 'both'). default: two-tailed
%                   if there is no a priori expectation, then use 'tail' = 'both'
%                   if the a priori expectation is that a > b, then use 'tail' = 'right'
%                   if the a priori expectation is that a < b, then use 'tail' = 'left'
%     dim:        (optional) dimension across which to test. default: 1
%
%   Output:
%     h:          significance (1 or 0)
%     pval:       p-value of permutation test. Discard H0 if pval is small.
%     diff:       mean(a)-mean(b)
%     diff_null:  the permuted null distribution
%
% RL van den Brink, 2019

%% check the input

if isempty(varargin) || nargin == 1
    b = zeros(size(a));
    npermutes = 1000;
    pthresh   = 0.05;
    tail = 'both';
    dim = 1;
else    
    if nargin == 0
        error('not enough input arguments')
    elseif nargin == 2
        npermutes = 1000;
        pthresh   = 0.05;
        tail = 'both';
        b = varargin{1};
        dim = 1;
    elseif nargin == 3
        pthresh   = 0.05;
        b = varargin{1};
        npermutes = varargin{2};
        tail = 'both';
        dim = 1;
    elseif nargin == 4
        b = varargin{1};
        npermutes = varargin{2};
        pthresh = varargin{3};
        tail = 'both';
        dim = 1;
    elseif nargin == 5
        b = varargin{1};
        npermutes = varargin{2};
        pthresh = varargin{3};
        tail = varargin{4};
        dim = 1;
    elseif nargin > 6
        error('too many input arguments')
    else
        b = varargin{1};
        npermutes = varargin{2};
        pthresh = varargin{3};
        tail = varargin{4};
        dim = varargin{5};
    end
end

%get rid of possible singleton dimensions
a = squeeze(a);
b = squeeze(b);

if (length(a) ~= length(b)) && length(b) == 1
    b = zeros(size(a)) + b;
elseif (length(a) ~= length(b)) && length(b) ~= 1
    error('The data in a paired sample test must be the same size.')
end

if isempty(b); b                 = zeros(size(a)); end
if isempty(npermutes); npermutes = 1000;   end
if isempty(pthresh); pthresh     = 0.05;   end
if isempty(tail); tail           = 'both'; end
if isempty(dim); dim             = 1; end
if sum(strcmpi(tail,{'both'; 'left'; 'right'})) ~= 1; error('Unrecognized tail option'); end
if ~isnumeric(pthresh); error('Provide a numeric p threshold'); end

rv = 0;
%make sure the input distributions are row vectors if necessary
if (length(size(a)) == 2 && length(size(b)) == 2) && (any(size(b) == 1) && any(size(a) == 1))
    rv = 1; %we're testing only vectors
    if size(a,2) > size(a,1); a = a'; end
    if size(b,2) > size(b,1); b = b'; end
end


%% test

%in case of vectors:
if rv %if we're only comparing the mean of vectors
    %compute difference in mean
    diff = mean(a)-mean(b);
    fprintf('Permuting:\t')
    %compute permuted null distribution of mean differences
    diff_null = zeros(npermutes,1);
    for permi = 1:npermutes       
        
        bnull = [a b];
        idx   = rand(size(a)) < 0.5;
        idx   = logical([idx 1-idx]);
        anull = bnull(idx);
        bnull(idx) = [];
        
        diff_null(permi) = mean(anull)-mean(bnull);
        if mod(permi,npermutes/10) == 0; fprintf([num2str((permi/npermutes)*100) '%%\t']); end
    end
    fprintf('Done!\n')
    %match the empirical matrix size to the null distribution
    diffm = ones(size(diff_null))*diff;
    
%if we're comparing the mean across a particular dimension of a matrix
else      
    s = size(a); %get datasize for data, later on
    s(dim) = [];
    diff = squeeze(mean(a,dim)-mean(b,dim));
    fprintf('Permuting:\t')
    %compute permuted null distribution of mean differences
    diff_null = zeros([npermutes prod(s)]);
    cnull = cat(dim,a,b); %concatinate across participants
    cnull = permute(cnull,[dim setdiff(1:length(size(cnull)),dim)]); %set participant dimension as the first
    for permi = 1:npermutes     
        
        idx  = rand(size(a,dim),1) < 0.5;
        idx  = logical([idx; idx-1]);
        anull = cnull;
        bnull = cnull;
        anull(~idx,:) = [];
        bnull(idx,:)  = [];
        d = mean(anull)-mean(bnull);
        
        diff_null(permi,:) = d(:);        
        if mod(permi,npermutes/10) == 0; fprintf([num2str((permi/npermutes)*100) '%%\t']); end
    end
    fprintf('Done!\n')
    diff_null = reshape(diff_null,[npermutes s]);
    if isvector(diff)
        if dim == 1
            diffm = repmat( diff, [npermutes,1] );
        else
            diffm = repmat( diff', [npermutes,1] );
        end
    else
        diffm = permute( repmat(diff,[ones(1,length(size(diff))) npermutes]), [length(size(diff_null)) 1:length(size(diff))]);
    end
end

%calculate p-value
if strcmpi(tail,'both')
    p = squeeze(1-sum(abs(diffm)>abs(diff_null))/npermutes);
elseif strcmpi(tail,'right')
    p = squeeze(1-sum(diffm>diff_null)/npermutes);
elseif strcmpi(tail,'left')
    p = squeeze(1-sum(diffm<diff_null)/npermutes);
end

%compare p to the alpha level
h = p < pthresh;

end