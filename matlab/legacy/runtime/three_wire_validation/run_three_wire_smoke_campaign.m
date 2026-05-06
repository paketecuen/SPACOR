function campaign = run_three_wire_smoke_campaign()
%RUN_THREE_WIRE_SMOKE_CAMPAIGN First 3-wire numerical validation campaign.

root = fileparts(mfilename('fullpath'));
addpath(root);
config = nv3_default_config();
if ~exist(config.resultsDir, 'dir')
    mkdir(config.resultsDir);
end

cases = {
    'delta_parallel_clean', @() nv3_synthetic_delta_parallel_case(config, struct('noiseSnrDb', Inf));
    'delta_parallel_snr60', @() nv3_synthetic_delta_parallel_case(config, struct('noiseSnrDb', 60));
    'wye_series_clean', @() nv3_synthetic_wye_series_case(config, struct('noiseSnrDb', Inf));
    'wye_series_snr60', @() nv3_synthetic_wye_series_case(config, struct('noiseSnrDb', 60));
    };

summaryRows = {};
parameterTables = cell(size(cases, 1), 1);
results = struct();

for idx = 1:size(cases, 1)
    caseName = cases{idx, 1};
    data = cases{idx, 2}();
    result = nv3_identify_three_wire_windowed(data, config, struct( ...
        'harmonics', 'adaptive', ...
        'windowCycles', config.windowCycles, ...
        'hopCycles', config.hopCycles));

    results.(caseName).data = data;
    results.(caseName).result = result;
    parameterTables{idx} = add_case_column(result.summary, caseName);

    userData = result.summary.Properties.UserData;
    maxRelError = max(result.summary.RelativeError(isfinite(result.summary.RelativeError)));
    if isempty(maxRelError)
        maxRelError = NaN;
    end
    summaryRows(end + 1, :) = {caseName, data.model, data.noiseSnrDb, ...
        userData.ValidCoverage, userData.MedianCondition, ...
        userData.MedianEquationResidual, userData.MedianPowerResidual, ...
        maxRelError}; %#ok<AGROW>

    writetable(result.windows, fullfile(config.resultsDir, [caseName '_windows.csv']));
    writetable(result.summary, fullfile(config.resultsDir, [caseName '_parameters.csv']));
    nv3_plot_parameter_windows(result, data, ...
        fullfile(config.resultsDir, [caseName '_parameters.png']), ...
        struct('truthMargin', config.plotTruthMargin));
end

summary = cell2table(summaryRows, 'VariableNames', {'Case', 'Model', 'NoiseSnrDb', ...
    'ValidCoverage', 'MedianCondition', 'MedianEquationResidual', ...
    'MedianPowerResidual', 'MaxRelativeError'});
parameters = vertcat(parameterTables{:});

writetable(summary, fullfile(config.resultsDir, 'three_wire_smoke_summary.csv'));
writetable(parameters, fullfile(config.resultsDir, 'three_wire_smoke_parameters.csv'));
save(fullfile(config.resultsDir, 'three_wire_smoke_campaign.mat'), ...
    'summary', 'parameters', 'results', 'config');

campaign = struct();
campaign.summary = summary;
campaign.parameters = parameters;
campaign.results = results;
campaign.config = config;

disp(summary);
disp(parameters);

end

function tableOut = add_case_column(tableIn, caseName)
tableOut = tableIn;
tableOut.Case = repmat({caseName}, height(tableOut), 1);
tableOut = movevars(tableOut, 'Case', 'Before', 1);
end
