#ifndef SIMULATOR_PILOT_H
#define SIMULATOR_PILOT_H

#include <cmath>
#include <complex>
#include <vector>
#include <utility>
#include <iostream>
#include <omp.h>
#include <Eigen/Dense>
#include "simulator_base.h"

class SimulatorPilot : public SimulatorBase {
public:
    SimulatorPilot(const SimulationParameters &params) : SimulatorBase(params) {}

    /**
     * 数値計算実験
     * @return パイロットシンボルのみをもちいたビット誤り率のシミュレーション値
     */
    double getBER_Simulation_only_pilot()
    {
        int totalErrorCount = 0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.equalizeByPilotAndDemodulate(transmitter_.getX());
            totalErrorCount += receiver_.getBitErrorCount(transmitter_.getTxData());
        }
        return (double)totalErrorCount / ((double)NUMBER_OF_TRIAL * (double)params_.NUMBER_OF_BIT * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

    /**
     * チャネル推定のMSEを計算するシミュレーション
     * @return パイロットシンボルのみを用いたMSEのシミュレーション値
     */
    double getMSE_simulation_only_pilot()
    {
        double totalSquaredError = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.equalizeByPilotAndDemodulate(transmitter_.getX());
            totalSquaredError += receiver_.getMSE();
        }
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

    double get_h_MSE_Simulation_during_pilot()
    {
        double totalSquaredError = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_initial_h(transmitter_.getX());
            totalSquaredError += receiver_.getMSE_during_pilot();
        }
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    double get_h_MSE_Simulation_during_pilot_RaghavendraGAIC()
    {
        double totalSquaredError = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_initial_h_RaghavendraGAIC(transmitter_.getX());
            totalSquaredError += receiver_.getMSE_during_pilot();
        }
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    double get_H_est_MSE_Simulation_during_pilot()
    {
        double totalSquaredError = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_pilot(transmitter_.getX());
            totalSquaredError += receiver_.getMSE_during_pilot();
        }
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * 真のパスモデルと既知の雑音分散を使った ML 推定の H のMSE を計算する
     */
    double get_H_MSE_with_known_model_and_noise_during_pilot()
    {
        double total_H_squared_error = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_known_model_and_noise(transmitter_.getX(), noiseSD_ * noiseSD_);
            total_H_squared_error += receiver_.getMSE_during_pilot();
        }
        return total_H_squared_error / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * パスモデルがわからないが16パスあると仮定してインパルス応答を推定する方法の H のMSE を計算する
     */
    double get_H_MSE_with_16paths_during_pilot()
    {
        double total_H_squared_error = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_16paths(transmitter_.getX());
            total_H_squared_error += receiver_.getMSE_during_pilot();
        }
        return total_H_squared_error / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * Raghavendraの論文で提案されているAICでモデル選択をした場合の推定
     */
    double get_H_MSE_with_RaghavendraAIC_during_pilot()
    {
        double total_H_squared_error = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_RaghavendraAIC(transmitter_.getX());
            total_H_squared_error += receiver_.getMSE_during_pilot();
        }
        return total_H_squared_error / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * 真のパスモデルと既知の雑音分散を使った ML 推定の H のMSE をフレーム先頭 (l=0) のみで計算する
     */
    double get_H_MSE_with_known_model_and_noise_at_l0()
    {
        double total_H_squared_error = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_known_model_and_noise(transmitter_.getX(), noiseSD_ * noiseSD_);
            total_H_squared_error += receiver_.getMSE_at_l0();
        }
        return total_H_squared_error / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * AICによるモデル選択の正答率を計算するシミュレーション
     */
    double getAICAccuracy_pilot()
    {
        int successCount = 0;
        #pragma omp parallel reduction(+:successCount)
        {
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 
            Channel local_channel(local_params, W_master_);
            Transmitter local_transmitter(local_params);
            Receiver local_receiver(local_params, W_master_);

            #pragma omp for
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transmitter.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_receiver.setY_(local_channel.getH(), local_transmitter.getX(), noiseSD_);
                local_receiver.est_H_by_initial_h(local_transmitter.getX()); 

                if (local_receiver.checkAICAccuracy(local_params.pathMask)) {
                    successCount++;
                }
            }
        }
        return (double)successCount / (double)NUMBER_OF_TRIAL;
    }

    /**
     * AICの評価指標（F値と正答率）を計算する
     */
    std::pair<double, double> getAIC_Metrics_pilot()
    {
        long long total_TP = 0; 
        long long total_FP = 0; 
        long long total_FN = 0; 
        long long total_TN = 0; 

        #pragma omp parallel reduction(+:total_TP, total_FP, total_FN, total_TN)
        {
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 
            Channel local_channel(local_params, W_master_);
            Transmitter local_transmitter(local_params);
            Receiver local_receiver(local_params, W_master_);

            #pragma omp for
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transmitter.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_receiver.setY_(local_channel.getH(), local_transmitter.getX(), noiseSD_);
                local_receiver.est_H_by_initial_h(local_transmitter.getX()); 

                std::vector<int> estMask = local_receiver.getEstimatedPathMask();
                
                for(int q=0; q<local_params.Q_; ++q) {
                    bool isTrue = (local_params.pathMask[q] == 1);
                    bool isEst  = (estMask[q] == 1);

                    if (isTrue && isEst)       total_TP++;
                    else if (!isTrue && isEst) total_FP++;
                    else if (isTrue && !isEst) total_FN++;
                    else                       total_TN++;
                }
            }
        }

        long long denominator_prec = total_TP + total_FP;
        double precision = (denominator_prec > 0) ? (double)total_TP / denominator_prec : 0.0;

        long long denominator_rec = total_TP + total_FN;
        double recall = (denominator_rec > 0) ? (double)total_TP / denominator_rec : 0.0;

        double numerator_f = 2.0 * precision * recall;
        double denominator_f = precision + recall;
        double f_measure = (denominator_f > 0) ? numerator_f / denominator_f : 0.0;

        long long total = total_TP + total_FP + total_FN + total_TN;
        double accuracy = (total > 0) ? (double)(total_TP + total_TN) / total : 0.0;

        return {f_measure, accuracy};
    }

    /**
     * 各試行でランダムなパス構成を生成し、平均MSEを計算する (Mode 40用)
     */
    double getMSE_RandomPath_Mode12_Simulation()
    {
        double totalSquaredError = 0.0;
        uniform_int_distribution<> dist;
        dist.init(0, 1, params_.seed);

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateRandomPathFrequencyResponse(fd_Ts_, dist);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_initial_h(transmitter_.getX());
            totalSquaredError += receiver_.getMSE_during_pilot();
        }
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * Mode 42: ランダムパスモデルによる平均MSE (Raghavendra AIC を全探索で適用)
     */
    double getMSE_RandomPath_RaghavendraAIC_Simulation()
    {
        double totalSquaredError = 0.0;
        uniform_int_distribution<> dist;
        dist.init(0, 1, params_.seed);

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateRandomPathFrequencyResponse(fd_Ts_, dist);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_RaghavendraAIC(transmitter_.getX());
            totalSquaredError += receiver_.getMSE_during_pilot();
        }
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * Mode 42: ランダムパスモデルによる平均MSE (Raghavendra AIC を全探索で適用) - 2
     */
    double getMSE_RandomPath_RaghavendraAIC_Simulation2()
    {
        double totalSquaredError = 0.0;
        uniform_int_distribution<> dist;
        dist.init(0, 1, params_.seed);

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateRandomPathFrequencyResponse(fd_Ts_, dist);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_RaghavendraAIC2(transmitter_.getX());
            totalSquaredError += receiver_.getMSE_during_pilot();
        }
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * ステップ 2: AIC 8パス総当たりによるモデル選択正答率のシミュレーション
     */
    double getExhaustiveAICAccuracy_8paths_Simulation()
    {
        int successCount = 0;
        uniform_int_distribution<> dist;
        dist.init(0, 1, params_.seed);

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            std::vector<int> true_mask_8(8);
            bool all_zero = true;
            for (int q = 0; q < params_.Q_; q++) {
                if (q < 8) {
                    int val = dist();
                    params_.pathMask[q] = val;
                    true_mask_8[q] = val;
                    if (val == 1) all_zero = false;
                } else {
                    params_.pathMask[q] = 0;
                }
            }

            if (all_zero) {
                int force_idx = tri % 8;
                params_.pathMask[force_idx] = 1;
                true_mask_8[force_idx] = 1;
            }

            channel_.updateProfile();
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);

            std::vector<int> est_mask_8 = receiver_.findBestMaskByExhaustiveAIC_8paths(transmitter_.getX());

            if (est_mask_8 == true_mask_8) {
                successCount++;
            }

            if (NUMBER_OF_TRIAL >= 10 && (tri + 1) % (NUMBER_OF_TRIAL / 10) == 0) {
                std::cout << "\rExhaustive AIC Simulation Progress: " << (int)((double)(tri + 1) / NUMBER_OF_TRIAL * 100.0) << "%" << std::flush;
            }
        }
        std::cout << std::endl;

        return (double)successCount / (double)NUMBER_OF_TRIAL;
    }

    /**
     * ステップ 2: AIC 8パス総当たりによるモデル選択後のインパルス応答推定MSEシミュレーション (固定パス用)
     */
    double getMSE_ExhaustiveAIC_8paths_fixedMask_Simulation()
    {
        double totalSquaredError = 0.0;
        channel_.updateProfile();

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);

            std::vector<int> best_mask_8 = receiver_.findBestMaskByExhaustiveAIC_8paths(transmitter_.getX());

            std::vector<int> full_mask(params_.Q_, 0);
            for (int q = 0; q < 8; ++q) full_mask[q] = best_mask_8[q];

            receiver_.est_H_by_given_mask(transmitter_.getX(), full_mask);
            totalSquaredError += receiver_.getMSE_during_pilot();

            if (NUMBER_OF_TRIAL >= 10 && (tri + 1) % (NUMBER_OF_TRIAL / 10) == 0) {
                std::cout << "\rExhaustive AIC MSE Progress: " << (int)((double)(tri + 1) / NUMBER_OF_TRIAL * 100.0) << "%" << std::flush;
            }
        }
        std::cout << std::endl;

        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * ステップ 2 (Raghavendra版): Raghavendra AIC 8パス総当たりによるモデル選択後のインパルス応答推定MSEシミュレーション (固定パス用)
     */
    double getMSE_ExhaustiveRaghavendraAIC_8paths_fixedMask_Simulation()
    {
        double totalSquaredError = 0.0;
        channel_.updateProfile();

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);

            std::vector<int> best_mask_8 = receiver_.findBestMaskByExhaustiveRaghavendraAIC_8paths(transmitter_.getX());

            std::vector<int> full_mask(params_.Q_, 0);
            for (int q = 0; q < 8; ++q) full_mask[q] = best_mask_8[q];

            receiver_.est_H_by_given_mask(transmitter_.getX(), full_mask);
            totalSquaredError += receiver_.getMSE_during_pilot();

            if (NUMBER_OF_TRIAL >= 10 && (tri + 1) % (NUMBER_OF_TRIAL / 10) == 0) {
                std::cout << "\rExhaustive Raghavendra AIC MSE Progress: " << (int)((double)(tri + 1) / NUMBER_OF_TRIAL * 100.0) << "%" << std::flush;
            }
        }
        std::cout << std::endl;

        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * Mode 46: ランダムパスモデルによる平均MSE (真のパスマスクが既知として推定)
     */
    double getMSE_RandomPath_KnownMask_Simulation()
    {
        double totalSquaredError = 0.0;
        uniform_int_distribution<> dist;
        dist.init(0, 1, params_.seed);

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateRandomPathFrequencyResponse(fd_Ts_, dist);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);

            std::vector<int> true_mask(params_.Q_, 0);
            const Eigen::VectorXd& xi = channel_.getXi();
            for (int q = 0; q < params_.Q_; q++) {
                if (xi(q) > 0.0) {
                    true_mask[q] = 1;
                }
            }

            receiver_.est_H_by_given_mask(transmitter_.getX(), true_mask);
            totalSquaredError += receiver_.getMSE_during_pilot();
        }
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    /**
     * [Mode 39] 16パス推定時の平均インパルス応答電力 (l=0) を横軸 q で出力
     */
    void saveAverageEstimatedImpulseResponseByQ_16paths(const std::string& filename)
    {
        std::vector<double> avg_power_q(params_.Q_, 0.0);
        std::cout << "Calculating average ESTIMATED impulse response (16 paths) over " << NUMBER_OF_TRIAL << " trials..." << std::endl;

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_16paths(transmitter_.getX());

            Eigen::VectorXcd h_est = receiver_.getEstimatedPathCoefficients();

            for (int q = 0; q < params_.Q_; q++)
            {
                double power = std::norm(h_est(q)); 
                avg_power_q[q] += power;
            }

            if (NUMBER_OF_TRIAL >= 10 && ((tri + 1) % (NUMBER_OF_TRIAL / 10) == 0)) {
                std::cout << "\rProgress: " << (int)((double)(tri + 1) / NUMBER_OF_TRIAL * 100.0) << "%" << std::flush;
            }
        }
        std::cout << std::endl;

        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }
        ofs << "q,AveragePower,AverageMagnitude,Power_dB" << std::endl;

        for (int q = 0; q < params_.Q_; q++)
        {
            double mean_power = avg_power_q[q] / (double)NUMBER_OF_TRIAL;
            double mean_magnitude = std::sqrt(mean_power);
            double power_db = 10.0 * std::log10(mean_power + 1e-20);

            ofs << q << "," << mean_power << "," << mean_magnitude << "," << power_db << std::endl;
        }
        ofs.close();
    }

    /**
     * @brief Mode 47: AIC vs Raghavendra GAIC の単一試行評価
     */
    std::pair<std::vector<double>, std::vector<double>> getAICvsQ_SingleTrial_Simulation(double fd_Ts, double EbN0dB)
    {
        setDopplerFrequency(fd_Ts);
        setNoiseSD(EbN0dB);

        uniform_int_distribution<> dist;
        dist.init(0, 1, params_.seed);

        transmitter_.setX_();
        channel_.generateRandomPathFrequencyResponse(fd_Ts_, dist);
        receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);

        return receiver_.calculateAICvsQ(transmitter_.getX());
    }

    /**
     * @brief Mode 48: AIC vs Raghavendra GAIC の複数回平均評価
     */
    std::pair<std::vector<double>, std::vector<double>> getAICvsQ_Average_Simulation(double fd_Ts, double EbN0dB)
    {
        setDopplerFrequency(fd_Ts);
        setNoiseSD(EbN0dB);

        std::vector<double> sum_aic(params_.Q_, 0.0);
        std::vector<double> sum_gaic(params_.Q_, 0.0);

        uniform_int_distribution<> dist;
        dist.init(0, 1, params_.seed);

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateRandomPathFrequencyResponse(fd_Ts_, dist);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);

            auto [aic, gaic] = receiver_.calculateAICvsQ(transmitter_.getX());
            for (int q = 0; q < params_.Q_; ++q) {
                sum_aic[q] += aic[q];
                sum_gaic[q] += gaic[q];
            }

            if (NUMBER_OF_TRIAL >= 10 && (tri + 1) % (NUMBER_OF_TRIAL / 10) == 0) {
                std::cout << "\rAIC vs GAIC Simulation Progress: " << (int)((double)(tri + 1) / NUMBER_OF_TRIAL * 100.0) << "%" << std::flush;
            }
        }
        std::cout << std::endl;

        std::vector<double> avg_aic(params_.Q_);
        std::vector<double> avg_gaic(params_.Q_);
        for (int q = 0; q < params_.Q_; ++q) {
            avg_aic[q] = sum_aic[q] / (double)NUMBER_OF_TRIAL;
            avg_gaic[q] = sum_gaic[q] / (double)NUMBER_OF_TRIAL;
        }

        return {avg_aic, avg_gaic};
    }
};

#endif // SIMULATOR_PILOT_H
