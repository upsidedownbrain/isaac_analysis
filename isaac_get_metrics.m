function [isaac_struct] = isaac_get_metrics(tseries, rois, varargin)
% isaac_struct = isaac_get_metrics(tseries, rois, 'name1', value1,...)
% extracts the isaac metrics. The inputs have to be already formatted in a
% specific way (see description below), which can be done for nifti files
% with isaac_adapt_files()
%
% inputs:
%     (mandatory)
%
%     tseries:
%         Data matrix with the fMRI data, in column vector form (i.e. each
%         column is a time-series). This is easy to do, but a nifti file
%         can be converted to this format using isaac_adapt_files()
%         function.
%
%     rois:
%         Structure containing the regions to use, with the following
%         fields:
%           'names' (optional)
%              Cell array with the names of all the regions.
%           'idx'
%              Cell array with the indices of all the regions. Indices must
%              be integers (1-based) pointing to columns of tseries.
%           'as_x' and 'as_y'
%              Regions to be interpreted as the 'X' region and the 'Y'
%              region for shared metrics. They must be a vector of integers
%              pointing to cells in idx and names. If none are defined, an
%              'all for all' analysis will be run, i.e. all regions will be
%              considered as 'x' and 'y'.
%
%
%     (optional in name-value pairs)
%
%     'do_fisher':
%         If do_fisher=true, averages of correlations (DCorr and
%         homogeneity) are taken in Fisher-transformed data:
%           corrs. -> [ r-to-z -> average -> z-to-r ] -> avg. corr.
%         If do_fisher is set to false (the default), then the correlations
%         are averaged directly. Fisher transform stabilizes sample
%         variance and reduces the bias when averaging correlations.
%
%     'bx_estimation_method':
%         It can be either 'var_ratio' (default) or 'midpoint'. In method
%         'midpoint' Bx is the midpoint between its minimum and maximum
%         possible values (i.e. the validity interval that ensures that the
%         model is consistent). In method 'var_ratio' Bx is estimated so
%         that the variance of the shared part in each region is in the
%         same ratio than the homogeneous variance. t can be shown that
%         this value is always within the validity interval and empirically
%         appears more accurate.
% 
%      'do_force_hom_nonnegative':
%          Crop negative values of Hom and HVar to zero (by default false).
%          Since some homogeneous activity is in general present in fMRI
%          data, the expected values of HVar and Hom of any region are
%          non-negative. In some cases, however, they might be negative due
%          to sampling variability if (i) the homogeneous activity is too
%          low compared to the unstructured activity (ii) there are too few
%          voxels in the region. This leads to uninterpretable results
%          (negative estimated variance) of the inferential metrics at the
%          subject level, but unbiased estimators if averaging at a group
%          level. If do_force_hom_nonnegative = 1, negative values are set
%          to zero, which implies that homogeneous activity is not detected
%          (HVar=0, UVar=Var), and neither is shared/independent activity
%          with other regions (SVar = IVar = 0). This is meaningful at a
%          subject level, but might produce a positive bias when averaging
%          at a subject level.
%
%     'do_compute_full_matrix':
%          Debug only. By default do_compute_full_matrix=false. This option
%          is only for analyses in which some regions considered as 'x' are
%          also considered as 'y' (e.g. all-versus-all). For debugging one
%          can compute shared metrics between all x-y pairs
%          (do_compute_full_matrix=true), while in general repeated pairs
%          are omited (e.g. DCorr(ROI_A, ROI_A)).
% 
%     'do_get_full_matrix':
%          This option is for analyses in which some regions considered as
%          'x' are also considered as 'y' (e.g. all-versus-all) and when
%          do_compute_full_matrix is false. By default (do_get_full_matrix
%          = true), the not computed values are "faked" as if they would
%          have been computed (e.g. DCorr(ROI_A, ROI_A) = 1, and
%          DCorr(ROI_A, ROI_B) = DCorr(ROI_B, ROI_A) ). Otherwise the
%          not-computed values are 0 in the output matrix.
%          
%     'file_out':
%         If speficied, the metrics will be saved as a .mat file to the
%         file path speficied by file_out.
%
%
%
% outputs:
%         isaac_struct
%             struct with all the computed metrics. The field 'descriptive'
%             contains the descriptive metrics, i.e. homogeneity, variance
%             and distant correlation. The field 'inferential' contains the
%             variances of each component in the ISAAC model, i.e.
%             homogeneous, unstructured variances, and shared and
%             independent variances, as well as the estimated balance
%             coefficient Bx. The field 'info' contains useful information
%             such as the number of regions, the number of voxels from each
%             region and the number of time points that were used to
%             extract the timepoints.
%
%             The local metrics are shaped as a 1-D vector, separately for
%             'x' and 'y' regions. The shared metrics are shaped as a 2-D
%             matrix, in which rows correspond to the 'X' region and the
%             row to the 'Y' region. This is important to interpret the Bx
%             coefficient, as well as the independent variance, for
%             example.
%
%             Additionally, the field 'raw' contains the average values
%             from the matrix, and the field 'meansignals' contains
%             variance, covariance and correlation from mean signals.
% 
% 
% 
% Notes:
%
%   The concept of homogeneity and homogeneous variance (HVar) is not
% well defined for one-voxel regions. In this special case, the region will
% be considered to have purely homogeneous activity. This implies, Hom = 1,
% HVar = Var, UVar = 0.
%
%   The ISAAC model assumes that there is some homogeneous activity
% within a region. If there is absolutely no homogeneous signal (i.e. Hom
% and HVar are theoretically zero) in a region there is a chance that
% Homogeneity and HVar estimates have negative values, which doesn't make
% sense. This is virtually impossible in real data, as some smoothness is
% always expected, but if the region is very small and homogeneous variance
% is very low related to noise (of the order of SNR<0.1), it might happen.


% --------- parse arguments -----------------------------------------------
if rem(length(varargin), 2)
    warning('Odd number of optional arguments. Optional arguments must be in name-value pairs');
end

% default values:
do_fisher    = false;
file_out   = '';
bx_estimation_method    = 'var_ratio';
do_get_full_matrix = true;
do_compute_full_matrix = false;
do_force_hom_nonnegative = false;

% If any values were specified, change the value:
for k = 1:2:length(varargin)
    name_ = varargin{k};
    switch name_
        case 'do_fisher';
            do_fisher = varargin{k+1};
        case 'do_get_full_matrix';
            do_get_full_matrix = varargin{k+1};
        case 'do_compute_full_matrix';
            do_compute_full_matrix = varargin{k+1};
        case 'bx_estimation_method'
            if ismember(varargin{k+1}, {'midpoint', 'var_ratio'})
                bx_estimation_method = varargin{k+1};
            else
                warning('bx_estimation_method %s not recognized, must be either {''midpoint'', ''var_ratio''}. ''%s'' value will be used', ...
                    varargin{k+1}, bx_estimation_method);
            end
        case 'do_force_hom_nonnegative';
            do_force_hom_nonnegative = varargin{k+1};
        case 'file_out';
            file_out = varargin{k+1};
        otherwise
            warning('Parameter ''%s'' not recognized, it will be ignored', name_);
    end
end
do_save = ~isempty(file_out);
% do_compute_full_matrix = true overrides do_get_full_matrix
do_get_full_matrix = do_get_full_matrix && (~do_compute_full_matrix);



% ------------- manage rois -----------------------------------------------
% if as_x and as_y not specified, do all versus all analysis
if ~isfield(rois, 'as_x') && ~isfield(rois, 'as_y')
    rois.as_x = 1:numel(rois.idx);
    rois.as_y = 1:numel(rois.idx);
end

% if the names are not specified, set some
if ~isfield(rois, 'names')
    rois.names = arrayfun(@(num) sprintf('roi%03.0f', num), 1:length(rois.idx), 'uni', false);
end



% ----------- compute the metrics -----------------------------------------
% Ensure tseries are demeaned:
tseries = tseries - repmat(mean(tseries, 1), size(tseries, 1), 1);

% basic data curation:
% remove zero-variation indices (they break dcorr if do_fisher==true)
idx_ok = cellfun(@(idx) idx(var(tseries(:,idx))>0), rois.idx, 'uni', false);
if any(cellfun(@(i1, i2) numel(i1)~=numel(i2), idx_ok, rois.idx))
    warning('Some regions contained voxels with no activity in tseries, these voxels were discarded.');
end
rois.idx = idx_ok;


isaac_struct = tseries_to_metrics(tseries,  rois, do_fisher, bx_estimation_method, do_get_full_matrix, do_compute_full_matrix, do_force_hom_nonnegative);

% add flags to metadata:
isaac_struct.info.bx_estimation_method = {bx_estimation_method};
isaac_struct.info.do_fisher_average_corrs = do_fisher;
isaac_struct.info.do_force_hom_nonnegative = do_force_hom_nonnegative;
% operations for statistics were not performed:
isaac_struct.info.fisher_transformed_corrs = false;
isaac_struct.info.log10_transformed_vars = false;




% ------------- save in mat file ------------------------------------------
% First make sure the file name is correct.
if do_save
    % if prefix is a dir, add file name. if it doesn't end with .mat, add
    % .mat.
    is_dir_file_out = (file_out(end)=='/') || (file_out(end)=='\') || isdir(file_out);
    is_mat_file_out = ~isempty(regexp(file_out, '\.mat$'));
    if is_dir_file_out
        dir_out = file_out;
        file_out = fullfile(file_out, sprintf('isaac_metrics_%s.mat', datestr(now, 'yyyymmdd_HHMMSSFFF')));
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
    save(file_out, '-struct', 'isaac_struct');
end



%##########################################################################
% Auxiliary functions
%##########################################################################

function [ isaac_struct ] = tseries_to_metrics(tseries, rois, do_fisher, bx_estimation_method, do_get_full_matrix, do_compute_full_matrix, do_force_hom_nonnegative)


if do_get_full_matrix
    % sort regions to be able to reconstruct the matrix later:
    [rois.as_x, is_rois_as_x] = sort(rois.as_x);
    [rois.as_y, is_rois_as_y] = sort(rois.as_y);
    % indices to sort back to user-defined order (a_unsorted =
    % a_sorted(ius))
    ius_rois_as_x = arrayfun(@(a) find(is_rois_as_x==a), 1:length(rois.as_x));
    ius_rois_as_y = arrayfun(@(a) find(is_rois_as_y==a), 1:length(rois.as_y));
end



% ------------ gather information -----------------------------------------
n_vols = size(tseries, 1);

info = struct;
% Indices from the image corresponding to the rois
info.idx   = rois.idx;
info.idx_x = rois.idx(rois.as_x);
info.idx_y = rois.idx(rois.as_y);

% Get "metadata" of the analysis
info.n_vols      = n_vols;
info.roi_names   = rois.names;
info.roi_names_x = rois.names(rois.as_x);
info.roi_names_y = rois.names(rois.as_y);
info.n_voxels    = cellfun(@numel, rois.idx);
info.n_voxels_x  = cellfun(@numel, rois.idx(rois.as_x));
info.n_voxels_y  = cellfun(@numel, rois.idx(rois.as_y));

n_rois = numel(info.idx);
n_rois_x = numel(rois.as_x);
n_rois_y = numel(rois.as_y);




% -------------------------------------------------------------------------
% First compute the "raw" metrics, i.e. the average of correlation and
% covariance matrices at a voxel level, and teh corrs/covs/variances of
% mean signals.
% 
% Then use raw_to_isaac() to estimate and arrange the inferential and
% descriptive metrics in the final output structure.
% -------------------------------------------------------------------------




% -------------------- initialize  metrics --------------------------------
% raw metrics
raw = struct;
[s_mean, s, c, r] = deal(zeros(n_rois,1));
% [raw.s_mean_x, raw.s_x, raw.c_xx, raw.r_xx] = deal(nan(n_rois_x,1));
% [raw.s_mean_y, raw.s_y, raw.c_yy, raw.r_yy] = deal(nan(n_rois_y,1));
[raw.c_xy, raw.r_xy, raw.r_mean_xy] = deal(zeros(n_rois_x, n_rois_y));
% degs of freedom for covs:
df_cov = (info.n_vols-1);



% -------- Compute local metrics ------------------------------------------
for k_x = 1:n_rois
    idx = info.idx{k_x};
    
    % extract submatrix (get corr from cov, don't compute twice)
    cov_submatrix = (tseries(:,idx)'*tseries(:,idx)) / df_cov;
    cor_submatrix = cov_submatrix./sqrt((diag(cov_submatrix)*diag(cov_submatrix)'));
    
    % For off-diagonal averaging. corr<1 is technically not the off-diag,
    % but we also remove (presumably artifactual) correlations = 1 (with
    % fisher r-to-z, only one element in cor_submatrix would drive the
    % average to 1)
    off_diag = cor_submatrix<1;
    % off_diag = ~eye(length(idx));
    
    % Average parts
    s_mean(k_x) = mean2(cov_submatrix);
    s(k_x)      = mean(diag(cov_submatrix));
    c(k_x)      = mean(cov_submatrix(off_diag));
    
    % correlation:
    if do_fisher
        r(k_x)     = tanh(mean(atanh(cor_submatrix(off_diag))));
    else
        r(k_x)     = mean(cor_submatrix(off_diag));
    end
end
% Special case for 1-voxel regions. These don't have off-diagonal elements.
% One voxel regions should have Hom=1, and HVar equal to the voxel's
% variance. For that we force local r to to 1, and HVar to the local
% variance.
c(info.n_voxels==1) = s(info.n_voxels==1);
r(info.n_voxels==1) = 1;

% detect c < 0, report
is_meancov_neg = c < 0;
is_meancor_neg = r < 0;
if sum(is_meancov_neg)>0 || sum(is_meancor_neg)>0
    warning_msg = ['Some regions have negative average within-region covariance and/or correlation. '...
        'This might render a negative homogeneity and non-interpretable inferential metrics (HVar, UVar, SVar, IVar). ' ...
        'The variance of homogeneous activity is too low to be detected or there are too few voxels'];
    if do_force_hom_nonnegative
        c(c<0) = 0;
        r(r<0) = 0;
        warning_msg = [warning_msg ...
            '\n do_force_hom_nonnegative is set to true, so negative values will be set to zero.'];
    end
    warning(sprintf(warning_msg));
end
% assign to x or y:
raw.s_mean_x = s_mean(rois.as_x);
raw.s_x      = s(rois.as_x);
raw.c_xx     = c(rois.as_x);
raw.r_xx     = r(rois.as_x);

raw.s_mean_y = s_mean(rois.as_y);
raw.s_y      = s(rois.as_y);
raw.c_yy     = c(rois.as_y);
raw.r_yy     = r(rois.as_y);


% -------- Compute shared metrics (x vs y) --------------------------------

% This saves time. Avoid nested loop iterating directly through all
% possible x and y combinations. For speed, remove repeated rois (A,A) and
% repeated pairs (A,B ~ B,A)
%  1. Take all possible roi pairs, and the indices:
[all_x, all_y]     = ndgrid(rois.as_x, rois.as_y);
[all_i_x, all_i_y] = ndgrid(1:n_rois_x, 1:n_rois_y);

% save the pairs that are equal, for later
is_rep_mat = (all_x == all_y);
all_xy = [all_x(:), all_y(:)];
all_i_xy = [all_i_x(:), all_i_y(:)];

if ~do_compute_full_matrix
    % 2. remove repetitions (A vs A)
    all_i_xy = all_i_xy(all_xy(:,2)~=all_xy(:,1),:);
    all_xy = all_xy(all_xy(:,2)~=all_xy(:,1),:);
    % 3. remove repeated pairs (B vs A if A vs B is computed)
    [all_xy, ia, ic] = unique(sort(all_xy, 2), 'rows', 'stable');
    all_i_xy = all_i_xy(ia,:);
end

for k_xy = 1:size(all_i_xy,1)
    k_x = all_i_xy(k_xy, 1);
    k_y = all_i_xy(k_xy, 2);
    idx_x = info.idx_x{k_x};
    idx_y = info.idx_y{k_y};
    
    % get cov and cor matrices:
    cov_submatrix = (tseries(:,idx_x)'*tseries(:,idx_y)) / df_cov;
    % same result, faster than corr(tseries(:,i_x), tseries(:,i_y));
    cor_submatrix = cov_submatrix./(std(tseries(:,idx_x))'*std(tseries(:,idx_y)));
    
    % for ISAAC metrics:
    raw.c_xy(k_x,k_y)  = mean2(cov_submatrix);
    
    if do_fisher
        raw.r_xy(k_x,k_y) = tanh(mean2(atanh(cor_submatrix)));
    else
        raw.r_xy(k_x,k_y) = mean2(cor_submatrix);
    end
    
    % Compute correlation of mean signal as cov/stds
    raw.r_mean_xy(k_x, k_y) = raw.c_xy(k_x, k_y) / sqrt(raw.s_mean_x(k_x)*raw.s_mean_y(k_y));
end



if do_get_full_matrix
    % if asked to, fill the not-computed elements of the matrices with the
    % same values that would have been computed if repetitions (A vs A) and
    % repeated pairs (B vs A and A vs B) would have been computed.
    
    % Duplicated pairs:
    rep_x = ismember(rois.as_x, rois.as_y);
    rep_y = ismember(rois.as_y, rois.as_x);
    raw.c_xy(rep_x, rep_y) = raw.c_xy(rep_x, rep_y) + raw.c_xy(rep_x, rep_y)';
    raw.r_xy(rep_x, rep_y) = raw.r_xy(rep_x, rep_y) + raw.r_xy(rep_x, rep_y)';
    raw.r_mean_xy(rep_x, rep_y) = raw.r_mean_xy(rep_x, rep_y) + raw.r_mean_xy(rep_x, rep_y)';
    
    % Repeated pairs (corr. is 1, cov w/ itself is the variance)
    raw.c_xy(is_rep_mat) = raw.s_mean_x(sum(is_rep_mat,2)>0);
    raw.r_xy(is_rep_mat) = 1;
    raw.r_mean_xy(is_rep_mat) = 1;
    
    % "unsort" back to user-defined order:
    raw.c_xy = raw.c_xy(ius_rois_as_x, ius_rois_as_y);
    raw.r_xy = raw.r_xy(ius_rois_as_x, ius_rois_as_y);
    raw.r_mean_xy = raw.r_mean_xy(ius_rois_as_x, ius_rois_as_y);
    raw.s_x = raw.s_x(ius_rois_as_x);
    raw.c_xx = raw.c_xx(ius_rois_as_x);
    raw.r_xx = raw.r_xx(ius_rois_as_x);
    raw.s_mean_x = raw.s_mean_x(ius_rois_as_x);
    raw.s_y = raw.s_y(ius_rois_as_y);
    raw.c_yy = raw.c_yy(ius_rois_as_y);
    raw.r_yy = raw.r_yy(ius_rois_as_y);
    raw.s_mean_y = raw.s_mean_y(ius_rois_as_y);
    
    rois.as_x = rois.as_x(ius_rois_as_x);
    rois.as_y = rois.as_y(ius_rois_as_y);
    
    info.roi_names_x = rois.names(rois.as_x);
    info.roi_names_y = rois.names(rois.as_y);
    info.n_voxels    = cellfun(@numel, rois.idx);
    info.n_voxels_x  = cellfun(@numel, rois.idx(rois.as_x));
    info.n_voxels_y  = cellfun(@numel, rois.idx(rois.as_y));
end



% Package the metrics and estimate inferential
isaac_struct = raw_to_isaac( raw, bx_estimation_method );

% Keep track of metadata to understand the dimensions
isaac_struct.info.n_time_points = info.n_vols;
isaac_struct.info.roi_names_x   = info.roi_names_x;
isaac_struct.info.roi_names_y   = info.roi_names_y;
isaac_struct.info.n_voxels_x    = info.n_voxels_x;
isaac_struct.info.n_voxels_y    = info.n_voxels_y;



function [ isaac_struct ] = raw_to_isaac( raw, bx_estimation_method )

% --------- store the "raw" averages from the matrices --------------------
isaac_struct.raw  = raw;

% --------- metrics extracted from the mean signals -----------------------
isaac_struct.meansignals = struct;
isaac_struct.meansignals.corr  = raw.r_mean_xy;
isaac_struct.meansignals.cov   = raw.c_xy;
isaac_struct.meansignals.var_x = raw.s_mean_x;
isaac_struct.meansignals.var_y = raw.s_mean_y;


% ------------  descriptive metrics ---------------------------------------
isaac_struct.descriptive = struct;
isaac_struct.descriptive.dcorr = raw.r_xy;
isaac_struct.descriptive.var_x = raw.s_x;
isaac_struct.descriptive.var_y = raw.s_y;
isaac_struct.descriptive.hom_x = raw.r_xx;
isaac_struct.descriptive.hom_y = raw.r_yy;


% -------------- inferential metrics --------------------------------------
isaac_struct.inferential = struct;

% homogeneous/unstructured decomposition (local)
isaac_struct.inferential.hvar_x   = raw.c_xx;
isaac_struct.inferential.hvar_y   = raw.c_yy;
isaac_struct.inferential.uvar_x   = raw.s_x - raw.c_xx;
isaac_struct.inferential.uvar_y   = raw.s_y - raw.c_yy;

% shared/independent decomposition:
[Bx, By] = estimate_bx( raw, bx_estimation_method );
[sz_x, sz_y] = size(raw.c_xy);

isaac_struct.inferential.svar_xy    = abs(raw.c_xy)./(Bx.*(1-Bx));
isaac_struct.inferential.svar_xyx   = abs(raw.c_xy).*(Bx./(1-Bx));
isaac_struct.inferential.svar_xyy   = abs(raw.c_xy).*((1-Bx)./Bx);
isaac_struct.inferential.ivar_xyx   = repmat(raw.c_xx, 1, sz_y, 1)  - abs(raw.c_xy).*(Bx./(1-Bx));
isaac_struct.inferential.ivar_xyy   = repmat(raw.c_yy', sz_x, 1, 1) - abs(raw.c_xy).*((1-Bx)./Bx);
isaac_struct.inferential.bx         = Bx;
isaac_struct.inferential.by         = By;


% Fix infs and nans. This is necessary to produce meaningful values if
% do_force_hom_nonnegative=true. If HVar = 0, SVar = 0, IVar = 0. Bx and By
% are set to 0, but they make no sense.
is_zero_x = raw.c_xx == 0;
is_zero_y = raw.c_yy == 0;
isaac_struct.inferential.svar_xy(is_zero_x, :) = 0;
isaac_struct.inferential.svar_xyx(is_zero_x, :) = 0;
isaac_struct.inferential.svar_xyy(is_zero_x, :) = 0;
isaac_struct.inferential.ivar_xyx(is_zero_x, :) = 0;
isaac_struct.inferential.ivar_xyy(is_zero_x, :) = 0;
isaac_struct.inferential.bx(is_zero_x, :) = 0;
isaac_struct.inferential.by(is_zero_x, :) = 0;


isaac_struct.inferential.svar_xy(:, is_zero_y) = 0;
isaac_struct.inferential.svar_xyx(:, is_zero_y) = 0;
isaac_struct.inferential.svar_xyy(:, is_zero_y) = 0;
isaac_struct.inferential.ivar_xyx(:, is_zero_y) = 0;
isaac_struct.inferential.ivar_xyy(:, is_zero_y) = 0;
isaac_struct.inferential.bx(:, is_zero_y) = 0;
isaac_struct.inferential.by(:, is_zero_y) = 0;

function [Bx, By] = estimate_bx( raw, bx_estimation_method )
switch bx_estimation_method
    case 'midpoint'
        % ----------------- midpoint of validity interval -----------------
        
        c   = abs(raw.c_xy);
        hx  =  repmat(abs(raw.c_xx(:)), 1, size(c,2));
        hy  =  repmat(abs(raw.c_yy(:)'), size(c,1), 1);
        
        Bx_min = 1./(1 + hy./c);
        Bx_max = 1./(1 + c./hx);
        Bx = .5*(Bx_min + Bx_max);
    case 'var_ratio'
        % ----------------- keep variance ratio ---------------------------
        %  The shared component is present in both X and Y, as Bx*s_xy and
        %  (1-Bx)*s_xy. This option estimates Bx so that:
        %      var(s_xy_x)/var(s_xy_y) = (Bx/(1-Bx))^2 = HVar(X)/HVar(Y)
        %  Which seems reasonable, and empirically it appears to give
        %  better estimations. It can be shown that this value is always
        %  within the range: The inequations (Bx_hat < Bx_min) and (Bx_hat
        %  > Bx_max), lead to the conclusion that correlation is greater
        %  than one (c_xy/sqrt(c_xx*c_yy) > 1).
        
        s_ratio = sqrt(abs(raw.c_xx(:))*(1./abs(raw.c_yy(:))'));
        Bx = s_ratio./(1 + s_ratio);
        Bx(isnan(raw.c_xy)) = nan;
end

% In any case, the sign of By is the sign of the covariance (if it is 0, By
% should have a value)
s = sign(raw.c_xy);
s(s==0) = 1;
By = (1-Bx).*s;
