classdef VoxelMapsTest < matlab.unittest.TestCase
%VOXELMAPSTEST Checks the voxel-map statistics against MATLAB's own tests,
%and the pipeline on a small made-up dataset (no images or SPM needed).
%   From this project's folder:  setup_voxel_maps; results = runtests('tests')
%   Needs the Statistics and Machine Learning Toolbox (as the original code did).

    methods (TestClassSetup)
        function addPaths(tc)
            here = fileparts(mfilename('fullpath'));
            run(fullfile(here, '..', 'setup_voxel_maps.m'));
            tc.assumeTrue(exist('ttest2', 'file') == 2, 'Statistics and Machine Learning Toolbox not available');
        end
    end

    methods (Test)
        function testTtestMatchesTtest2(tc)
            [L, y] = VoxelMapsTest.toy(70, 40, 1);
            [t, p] = vlsm_ttest(L, y);
            [t2, p2] = vlsm_ttest(L, y, 'equal');
            for v = 1:size(L, 2)
                a = y(~L(:, v)); b = y(L(:, v));
                if numel(a) < 2 || numel(b) < 2, tc.verifyTrue(isnan(t(v))); continue; end
                [~, pw, ~, sw] = ttest2(a, b, 'VarType', 'unequal');
                [~, pe, ~, se] = ttest2(a, b);
                tc.verifyEqual([t(v) p(v)], [sw.tstat pw], 'AbsTol', 1e-9);
                tc.verifyEqual([t2(v) p2(v)], [se.tstat pe], 'AbsTol', 1e-9);
            end
        end

        function testRanksumMatchesRanksum(tc)
            [L, y] = VoxelMapsTest.toy(70, 40, 2);
            y = round(3 * y);                                % ties
            [z, p] = vlsm_ranksum(L, y);
            for v = 1:size(L, 2)
                a = y(~L(:, v)); b = y(L(:, v));
                if isempty(a) || isempty(b), continue; end
                [pr, ~, st] = ranksum(a, b, 'method', 'approximate');
                tc.verifyEqual([z(v) p(v)], [st.zval pr], 'AbsTol', 1e-9);
            end
        end

        function testBayesUnitInformation(tc)
            % closed form for a two-group design with the corrected exponent
            [L, y] = VoxelMapsTest.toy(50, 5, 3);
            b = vlsm_bayes(L, y, 'unit_information');
            n = numel(y);
            for v = 1:5
                X = [ones(n, 1), double(L(:, v))];
                e = y - X * (X \ y);
                R2 = 1 - sum(e .^ 2) / sum((y - mean(y)) .^ 2);
                tc.verifyEqual(b(v), ((n - 2) / 2) * log(1 + n) - ((n - 1) / 2) * log(1 + n * (1 - R2)), 'AbsTol', 1e-10);
            end
        end

        function testJzsMatchesIntegral(tc)
            % with k = 1 this is jzs_corbf (Wetzels & Wagenmakers 2012)
            for n = [20 100 257]
                for r = [0 0.1 0.3 0.6]
                    f = @(g) exp(((n - 2) / 2) * log1p(g) - ((n - 1) / 2) * log1p((1 - r ^ 2) * g) ...
                        - 1.5 * log(g) - n ./ (2 * g));
                    ref = log(sqrt(n / 2) / gamma(0.5) * integral(f, 0, Inf));
                    tc.verifyEqual(bf_regression(r ^ 2, n, 1, 'jzs'), ref, 'AbsTol', 1e-6);
                end
            end
        end

        function testFdrIgnoresUntestedVoxels(tc)
            p = [0.001 0.002 0.035 0.2 0.5 NaN NaN]';
            [thr, sig] = fdr_threshold(p, 0.05, 'bh');
            tc.verifyEqual(thr, 0.002);                     % 0.035 > 3 * 0.05 / 5; NaNs are not tests
            tc.verifyEqual(sig, logical([1 1 0 0 0 0 0])');
            [thr, sig] = fdr_threshold([0.2 0.5]', 0.05);
            tc.verifyEqual(thr, 0); tc.verifyFalse(any(sig));
        end

        function testJointPermutationKeepsPairs(tc)
            % the difference of a patient's two scores is unchanged by
            % permuting whole patients, so its null maximum is constant
            Y = randn(RandStream('mt19937ar', 'Seed', 4), 30, 2);
            [~, null_max] = permutation_threshold(@(Yp) Yp(:, 1) - Yp(:, 2), Y, 50, 0.05, 1, 0);
            tc.verifyEqual(null_max, repmat(max(abs(Y(:, 1) - Y(:, 2))), 50, 1), 'AbsTol', 1e-12);
        end

        function testNormaliseScores(tc)
            s = RandStream('mt19937ar', 'Seed', 5);
            c = randn(s, 60, 2);
            y = 2 * c(:, 1) + randn(s, 60, 1);
            y(3) = NaN;
            r = normalise_scores(y, c, 'residual');
            ok = ~isnan(r);
            tc.verifyEqual(nnz(~ok), 1);
            tc.verifyLessThan(abs(corr(r(ok), c(ok, 1))), 1e-10);   % covariate removed
            z = normalise_scores(y, c);
            tc.verifyEqual(sort(z(ok)), norminv((1:59)' / 60), 'AbsTol', 1e-12);
        end

        function testInteractionMatchesAnovan(tc)
            s = RandStream('mt19937ar', 'Seed', 6);
            n = 60;
            L = rand(s, n, 3) < 0.4;
            F = rand(s, n, 2) + 0.1;
            y = randn(s, n, 1) + 1.5 * L(:, 1) .* F(:, 1);
            p = vlsm_interaction(L, y, F);
            model = [1 1 0; 1 0 0; 0 1 0; 0 0 1];
            for v = 1:3
                pa = anovan(y, [double(L(:, v)), F], 'model', model, 'continuous', [2 3], 'display', 'off');
                tc.verifyEqual(p(v), pa(1), 'AbsTol', 1e-8);
            end
        end

        function testPipelineFindsPlantedEffect(tc)
            D = VoxelMapsTest.dataset();
            cfg = config_voxel_maps();
            cfg.verbose = false;
            cfg.vlsm.methods = {'ttest', 'ranksum', 'bayes', 'lpd'};
            cfg.vlsm.fdr_q = 0.05;
            cfg.vlsm.n_perm = 50;
            R = vlsm_analysis(D, cfg);
            tc.verifyEqual(numel(R.maps), 12);
            m = R.maps(strcmp({R.maps.method}, 'ttest') & strcmp({R.maps.target}, 'Hand'));
            planted = ismember(R.voxels, 1:20);
            tc.verifyGreaterThan(nnz(m.fdr_mask & planted), 10);         % found most planted voxels
            tc.verifyLessThan(nnz(m.fdr_mask & ~planted), 5);            % few false positives
            tc.verifyGreaterThan(mean(m.stat(planted)), 2);              % lesion -> lower hand score
            b = R.maps(strcmp({R.maps.method}, 'bayes') & strcmp({R.maps.target}, 'Hand'));
            tc.verifyGreaterThan(mean(b.stat(planted)), mean(b.stat(~planted)));
        end

        function testCovariateComparisonRuns(tc)
            D = VoxelMapsTest.dataset();
            cfg = config_voxel_maps();
            cfg.verbose = false;
            cfg.vlsm.methods = {'ttest'};
            cfg.vlsm.fdr_q = 0.05;
            C = compare_covariate_control(D, cfg);
            tc.verifyEqual(height(C.summary), 9);            % 3 sets x 3 targets
            tc.verifyEqual(height(C.overlap), 6);            % 2 comparisons x 3 targets
        end

        function testPreprocessLesion(tc)
            X = false(91, 10, 20);
            X(20, 5, [5 10 15]) = true;                     % drawn every 5 slices, left of midline
            Y = preprocess_lesion(X, struct('gap_fill_slices', 6, 'one_hemisphere', true));
            tc.verifyFalse(any(any(any(Y(1:45, :, :)))));   % folded onto the other side
            tc.verifyTrue(all(squeeze(Y(72, 5, 5:15))));    % 91 + 1 - 20 = 72, gaps filled
        end
    end

    methods (Static)
        function [L, y] = toy(n, v, seed)
            s = RandStream('mt19937ar', 'Seed', seed);
            L = rand(s, n, v) < linspace(0.02, 0.6, v);
            y = randn(s, n, 1) - 0.8 * L(:, 1);
        end

        function D = dataset()
            s = RandStream('mt19937ar', 'Seed', 7);
            n = 120; V = 300;
            D.ids = arrayfun(@(k) sprintf('p%d', k), (1:n)', 'UniformOutput', false);
            D.lesions = rand(s, n, V) < 0.3;
            region = rand(s, n, 1) < 0.3;                   % voxels 1-20: one region, lesioned together
            D.lesions(:, 1:20) = repmat(region, 1, 20);
            D.voxels = (1:V)';
            D.dims = [V 1 1];
            D.score_names = {'Hand', 'Finger'};
            D.covariate_names = {'Age', 'TimePost', 'LesionVolume'};
            D.covariates = [60 + 10 * randn(s, n, 1), 12 * rand(s, n, 1), sum(D.lesions, 2)];
            hand = 18 - 2.5 * region + randn(s, n, 1);
            finger = 17 + randn(s, n, 1);
            D.scores = [hand, finger];
        end
    end
end
