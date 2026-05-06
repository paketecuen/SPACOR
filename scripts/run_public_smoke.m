function results = run_public_smoke()
%RUN_PUBLIC_SMOKE Run the public reproducibility smoke tests.
%
% This entry point uses synthetic fixtures and synthetic degradation campaigns
% only.  It does not require measurement datasets.

startup_spacor();
results = run_regression_tests();
fprintf('Public SPACOR smoke test completed.\n');
end
