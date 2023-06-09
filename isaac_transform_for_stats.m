function [ isaac_struct_out ] = isaac_transform_for_stats( isaac_metrics_in, varargin )
% isaac_struct_out = isaac_transform_for_stats( isaac_metrics_in, 'name1', value1,...)
%
% Transform the metrics for statistical testing. Variances are transformed
% by the log10 and correlations by Fisher's r-to-z transform (covariance
% and "raw" averages are left unchanged) The transformation is not
% necessary to perform statistical testing on the metrics, but it is common
% practice to transform correlations and variances in this way
%
% inputs:
%     (mandatory)
%
%     isaac_metrics_in:
%         Struct or path to a file generated by isaac_concatenate_metrics.
%
%     (optional in name-value pairs)
%
%     'file_out':
%         If speficied, the transformed  metrics will be saved to the file
%         path file_out (.mat).
%
%     'do_fisher_corrs' and 'do_log10_vars':
%         Two independent flags for the fisher r-to-z transformation of
%         correlations and log10 transformation of variances. By default
%         both are true.
%
% outputs:
%      isaac_struct_out
%          The file with transformed metrics.
%



% --------- parse arguments -----------------------------------------------
if rem(length(varargin), 2)
    warning('Odd number of optional arguments. Optional arguments must be in name-value pairs');
end

% default values:
file_out   = '';
do_fisher_corrs = true;
do_log10_vars   = true;

% If any values were specified, change the value:
for k = 1:2:length(varargin)
    name_ = varargin{k};
    switch name_
        case 'file_out';
            file_out = varargin{k+1};
        case 'do_fisher_corrs';
            do_fisher_corrs = varargin{k+1};
        case 'do_log10_vars';
            do_log10_vars = varargin{k+1};
        otherwise
            warning('Parameter ''%s'' not recognized, it will be ignored', name_);
    end
end
do_save = ~isempty(file_out);



% -------- detect the format of the input ---------------------------------
% allow only file name (of .mat) or isaaac struct
switch class(isaac_metrics_in)
    case 'char'
        
        % detect the input format
        [p,f,e] = fileparts(isaac_metrics_in);
        switch e
            case '.mat'
                isaac_struct_out = load(isaac_metrics_in);
            otherwise
                error('File extension ''%s'' for input file not recognized, it should be .mat', e);
        end
    case 'struct'
        % check that all I need is there:        
        fns = fieldnames(isaac_metrics_in);
        ok_fns = {'info', 'meansignals', 'descriptive', 'inferential'};
        if ~all(ismember(ok_fns, fns))
            error('Input struct should have fieldnames:\n\t%s\nbut it has fieldnames:\n\t%s', ...
                strjoin(ok_fns, ', '), strjoin(fns, ', '));
        end
        
        isaac_struct_out = isaac_metrics_in;
        
    otherwise
        error('Input should be a struct or a file path, as the output by isaac_concatenate_files, but instead is %s',...
            class(isaac_metrics_in));
end


% -------------------- do the thing ---------------------------------------
% variances 
if do_log10_vars
    is_transf = isaac_struct_out.info.log10_transformed_vars(:);
    % meanSignals
    for fn = {'var_x' 'var_y'}
        isaac_struct_out.meansignals.(fn{1})(:,:,~is_transf) = log10(isaac_struct_out.meansignals.(fn{1})(:,:,~is_transf));
    end
    % descriptive
    for fn = {'var_x' 'var_y'}
        isaac_struct_out.descriptive.(fn{1})(:,:,~is_transf) = log10(isaac_struct_out.descriptive.(fn{1})(:,:,~is_transf));
    end
    % inferential
    for fn = {'hvar_x' 'hvar_y' 'uvar_x' 'uvar_y' 'svar_xy' 'svar_xyx' 'svar_xyy' 'ivar_xyx' 'ivar_xyy'}
        % if Var(h_x) is near 0, negative values might arise. This is
        % problematic when computing log10, so they are cropped to zero
        % without asking.
        if any(isaac_struct_out.inferential.(fn{1})(:)<0)
            isaac_struct_out.inferential.(fn{1})(isaac_struct_out.inferential.(fn{1})<0) = 0;
            warning('some values of %s are less than zero. They are set to zero before computing logarithm to avoid complex numbers.', fn{1});
        end
        isaac_struct_out.inferential.(fn{1})(:,:,~is_transf) = log10(isaac_struct_out.inferential.(fn{1})(:,:,~is_transf));
    end
end

% correlations 
if do_fisher_corrs
    is_transf = isaac_struct_out.info.fisher_transformed_corrs(:);
    % mean signals
    for fn = {'corr'}
        isaac_struct_out.meansignals.(fn{1})(:,:,~is_transf) = atanh(isaac_struct_out.meansignals.(fn{1})(:,:,~is_transf));
    end
    % descriptive
    for fn = {'dcorr' 'hom_x' 'hom_y'}
        isaac_struct_out.descriptive.(fn{1})(:,:,~is_transf) = atanh(isaac_struct_out.descriptive.(fn{1})(:,:,~is_transf));
    end
end
isaac_struct_out.info.fisher_transformed_corrs(:) = do_fisher_corrs;
isaac_struct_out.info.log10_transformed_vars(:)   = do_log10_vars;



% ------------ save -------------------------------------------------------
if do_save
    % if file_out is a dir, add file name. if it doesn't end with .mat, add
    % .mat.
    is_dir_file_out = (file_out(end)=='/') || (file_out(end)=='\') || isdir(file_out);
    is_mat_file_out = ~isempty(regexp(file_out, '\.mat$'));
    if is_dir_file_out
        dir_out = file_out;
        file_out = fullfile(dir_out, sprintf('isaac_metrics_trans_%s.mat', datestr(now, 'yyyymmdd_HHMMSSFFF')));
        warning('Output file name ''%s'' is a directory. Data will be saved to ''%s''',...
            dir_out, file_out);
    else
        dir_out = fileparts(file_out);
        if ~is_mat_file_out
            file_out = [file_out '.mat'];
        end
    end
    
    % if output directory doesn't exist, create it.
    if ~exist(dir_out, 'dir')
        mkdir(dir_out);
    end
    save(file_out, '-struct', 'isaac_struct_out');
end

