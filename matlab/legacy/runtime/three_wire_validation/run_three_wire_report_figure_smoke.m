function outputs = run_three_wire_report_figure_smoke()
%RUN_THREE_WIRE_REPORT_FIGURE_SMOKE Generate Phase 8 report figures.

config = nv3_default_config();
config.duration = 0.25;
config.windowCycles = 0.5;
config.hopCycles = 0.1;
config.conditionLimit = 1e8;
config.rankTolerance = 1e-8;
config.sigmaFloor = 1e-8;
config.modelResidualFloor = 2.5e-3;
config.energyResidualFloor = 2.5e-3;
config.energyResidualFactor = 5;
config.selectionRelativeTolerance = 0.10;
config.adaptiveHarmonics = true;
config.harmonicSelectionMode = 'residual';
config.harmonicCandidates = [1 3 5 7 9 11 13];
config.activeThreshold = 2e-2;

outDir = fullfile(config.resultsDir, 'figures', ...
    'three_wire_report_smoke');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

cases = {
    'delta_glc_noise60', @() nv3_synthetic_delta_parallel_case(config, ...
    struct('model', 'delta_parallel_glc', 'noiseSnrDb', 60, ...
    'randomSeed', 3301))
    'wye_rlc_noise60', @() nv3_synthetic_wye_series_case(config, ...
    struct('model', 'wye_series_rlc', 'noiseSnrDb', 60, ...
    'randomSeed', 3302))
    'mixed_delta_wye_noise60', @() nv3_synthetic_mixed_delta_wye_case( ...
    config, struct('noiseSnrDb', 60, 'randomSeed', 3303))
    };

rows = cell(size(cases, 1), 6);
for idx = 1:size(cases, 1)
    name = cases{idx, 1};
    data = cases{idx, 2}();
    data.label = strrep(name, '_', ' ');
    report = nv3_characterize_three_wire_load(data, config, ...
        struct('candidateMode', 'mixed_blackbox'));
    outFile = fullfile(outDir, [name '_report.png']);
    nv3_plot_characterization_report(report, data, outFile, ...
        struct('truthMargin', 0.15));
    rows(idx, :) = {name, report.recommendation.status, ...
        report.recommendation.model, report.recommendation.selectedCoverage, ...
        report.recommendation.medianPowerResidual, outFile};
end

outputs = cell2table(rows, 'VariableNames', { ...
    'CaseName', 'Status', 'Model', 'Coverage', 'MedianPowerResidual', ...
    'Figure'});
disp(outputs);
end
