#ifndef SIMULATOR_MISC_H
#define SIMULATOR_MISC_H

#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <Eigen/Dense>
#include "simulator_base.h"

class SimulatorMisc : public SimulatorBase {
public:
    SimulatorMisc(const SimulationParameters &params) : SimulatorBase(params) {}

    /**
     * @brief 全試行回数における伝送路の平均電力を計算するシミュレーション
     * @return double 平均電力
     */
    double getAveragePower_simulation()
    {
        double totalPower = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            channel_.generateFrequencyResponse(fd_Ts_);
            totalPower += channel_.getAveragePower();
        }
        return totalPower / (double)NUMBER_OF_TRIAL;
    }

    void saveChannelMagnitudeResponseToCSV(std::ofstream& ofs, double fd_Ts)
    {
        channel_.generateFrequencyResponse(fd_Ts);
        const auto& H = channel_.getH();
        
        int K = params_.K_;
        int L = params_.L_;

        ofs << "l,k,Magnitude" << std::endl;

        for (int l = 0; l < L; ++l)
        {
            for (int k = 0; k < K; ++k)
            {
                double magnitude = std::abs(H(l, k));
                ofs << l << "," << k << "," << magnitude << std::endl;
            }
        }
    }

    /**
     * Mode34: 瞬時信号対雑音電力比 γ(ΔH) の平均を計算するシミュレーション
     */
    double getInstantaneousSNR_Mode34_simulation()
    {
        double totalGamma = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_initial_h(transmitter_.getX());

            Eigen::RowVectorXcd delta = receiver_.computeDeltaHRow(0);
            Eigen::VectorXd gamma_vec = receiver_.computeGammaFromDeltaH(delta, 0, noiseSD_, 1.0, -1.0, true);
            double gamma_k0 = 0.0;
            if (gamma_vec.size() > 0) gamma_k0 = gamma_vec(0);

            totalGamma += gamma_k0;
        }
        return totalGamma / (double)NUMBER_OF_TRIAL;
    }

    /**
     * Mode35: 瞬時信号対雑音電力比 γ(ΔH) の平均を計算するシミュレーション
     */
    double getInstantaneousSNR_Mode35_simulation()
    {
        double totalGamma = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);

            Eigen::RowVectorXcd delta = receiver_.computeDeltaRowFromPhotoX(0, noiseSD_);
            Eigen::VectorXd gamma_vec = receiver_.computeGammaFromDeltaH(delta, 0, noiseSD_, 1.0, -1.0, true);
            double gamma_k0 = 0.0;
            if (gamma_vec.size() > 0) gamma_k0 = gamma_vec(0);

            totalGamma += gamma_k0;
        }
        return totalGamma / (double)NUMBER_OF_TRIAL;
    }

    /**
     * Mode36: 瞬時信号対雑音電力比 γ(ΔH) の平均を計算するシミュレーション
     */
    double getInstantaneousSNR_Mode36_simulation()
    {
        double totalGamma = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);

            Eigen::VectorXd gamma_vec = receiver_.computeGamma_78(0, noiseSD_, 1.0, -1.0, true);
            double gamma_k0 = 0.0;
            if (gamma_vec.size() > 0) gamma_k0 = gamma_vec(0);

            totalGamma += gamma_k0;
        }
        return totalGamma / (double)NUMBER_OF_TRIAL;
    }

    /**
     * [Mode 23] 送信信号波形の出力実行関数
     */
    void runExportTxWaveform(int target_k, const std::string& filename)
    {
        SimulationParameters local_params = params_;
        Transmitter local_transmitter(local_params);
        local_transmitter.setX_();
        local_transmitter.exportTxSymbolTrace(target_k, filename);
    }

    /**
     * [Mode 24] フェージングを受けた信号波形の出力実行関数
     */
    void runExportFadedWaveform(int target_k, const std::string& filename)
    {
        SimulationParameters local_params = params_;
        Channel local_channel(local_params, W_master_);
        Transmitter local_transmitter(local_params);

        local_transmitter.setX_();
        local_channel.generateFrequencyResponse(fd_Ts_);
        local_transmitter.exportFadedSymbolTrace(target_k, filename, local_channel.getH());
    }

    /**
     * [Mode 25] チャネル応答の絶対値波形の出力実行関数
     */
    void runExportChannelMagnitude(int target_k, const std::string& filename)
    {
        SimulationParameters local_params = params_;
        Channel local_channel(local_params, W_master_);
        Transmitter local_transmitter(local_params);

        local_channel.generateFrequencyResponse(fd_Ts_);
        local_transmitter.exportChannelMagnitudeTrace(target_k, filename, local_channel.getH());
    }

    /**
     * [Mode 26] 特定の時刻 l における周波数応答の絶対値 |H(k, l)| を横軸 k で出力
     */
    void saveFrequencyResponseByK(int target_l, const std::string& filename)
    {
        channel_.generateFrequencyResponse(fd_Ts_);
        const auto& H = channel_.getH();
        
        std::ofstream ofs(filename);
        ofs << "k,Magnitude" << std::endl;
        for (int k = 0; k < params_.K_; ++k) {
            ofs << k << "," << std::abs(H(target_l, k)) << std::endl;
        }
        ofs.close();
    }

    /**
     * [Mode 27] 特定の時刻 l における真のインパルス応答 |h(q, l)| を横軸 q で出力
     */
    void saveImpulseResponseByQ(int target_l, const std::string& filename)
    {
        channel_.generateFrequencyResponse(fd_Ts_);
        const auto& h_true = channel_.get_h(); 

        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        ofs << "q,Magnitude" << std::endl;

        for (int q = 0; q < params_.Q_; ++q) {
            double magnitude = std::abs(h_true(target_l, q));
            ofs << q << "," << magnitude << std::endl;
        }
        ofs.close();
    }

    /**
     * [Mode 27修正版] 複数試行の平均インパルス応答 |h(q)| を横軸 q で出力
     */
    void saveAverageImpulseResponseByQ(int target_l, const std::string& filename)
    {
        std::vector<double> avg_power_q(params_.Q_, 0.0);
        std::cout << "Calculating average impulse response over " << NUMBER_OF_TRIAL << " trials..." << std::endl;

        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            channel_.generateFrequencyResponse(fd_Ts_);
            const auto& h_true = channel_.get_h();

            for (int q = 0; q < params_.Q_; q++)
            {
                double power = std::norm(h_true(target_l, q)); 
                avg_power_q[q] += power;
            }

            if (NUMBER_OF_TRIAL >= 10 && ((tri + 1) % (NUMBER_OF_TRIAL / 10) == 0)) {
                std::cout << "\rProgress: " << (int)((double)(tri + 1) / NUMBER_OF_TRIAL * 100.0) << "%" << std::flush;
            }
        }
        std::cout << std::endl;

        std::ofstream ofs(filename);
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
     * [Mode 28] パイロット区間 (l=0) における推定インパルス応答の保存
     */
    void saveEstimatedImpulseResponseToCSV(std::ofstream& ofs, double fd_Ts) {
        transmitter_.setX_();
        channel_.generateFrequencyResponse(fd_Ts);
        receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
        receiver_.est_H_by_initial_h(transmitter_.getX());

        Eigen::VectorXcd h_est = receiver_.getEstimatedPathCoefficients();

        ofs << "Path_Index,Real,Imag,Magnitude" << std::endl;
        for (int q = 0; q < h_est.size(); ++q) {
            double mag = std::abs(h_est(q));
            ofs << q << "," << h_est(q).real() << "," << h_est(q).imag() << "," << mag << std::endl;
        }
    }

    /**
     * ベイジアン・クラメル・ラオ下界 (BCRLB) の計算
     */
    double getTheoreticalCRLB_MSE()
    {
        double noiseVar = noiseSD_ * noiseSD_;
        const Eigen::VectorXd& xi = channel_.getXi();

        double crlb_h = 0.0;
        double sum_X_sq = (double)params_.K_; 

        for (int q = 0; q < params_.Q_; q++)
        {
            if (params_.pathMask[q] == 1)
            {
                double fisher_info = (sum_X_sq / noiseVar) + (1.0 / xi(q));
                crlb_h += 1.0 / fisher_info;
            }
        }
        return crlb_h;
    }

    /**
     * 誤差伝播則を用いた 周波数応答 (H) の CRLB (MSE理論下界) の計算
     */
    double getTheoreticalCRLB_H_MSE()
    {
        double noiseVar = noiseSD_ * noiseSD_;
        double sum_X_sq = (double)params_.K_; 

        Eigen::MatrixXcd C_h = Eigen::MatrixXcd::Zero(params_.Q_, params_.Q_);
        for (int q = 0; q < params_.Q_; q++)
        {
            if (params_.pathMask[q] == 1) 
            {
                C_h(q, q) = noiseVar / sum_X_sq; 
            }
        }

        Eigen::MatrixXcd C_H = W_master_ * C_h * W_master_.adjoint();
        double total_mse_H = C_H.trace().real();

        return total_mse_H / (double)params_.K_;
    }

    /**
     * 画像の最終式に基づく 周波数応答 (H) の CRLB (理論下界) の計算
     */
    double getTheoreticalCRLB_H_MSE_FinalForm()
    {
        double noiseVar = noiseSD_ * noiseSD_;
        double sum_X_sq = (double)params_.K_ * (double)params_.NUMBER_OF_PILOT; 

        int active_paths = 0;
        for (int q = 0; q < params_.Q_; q++)
        {
            if (params_.pathMask[q] == 1) 
            {
                active_paths++;
            }
        }
        return (double)active_paths * noiseVar / sum_X_sq; 
    }

    /**
     * インパルス応答（h）のMSEを計算するシミュレーション (パイロット区間)
     */
    double getImpulseResponseMSE_simulation()
    {
        double totalSquaredError = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.est_H_by_initial_h(transmitter_.getX()); 
            
            const Eigen::MatrixXcd& h_true_matrix = channel_.get_h();
            Eigen::VectorXcd h_true_l0 = h_true_matrix.row(0).transpose();
            Eigen::VectorXcd h_true_l1 = h_true_matrix.row(1).transpose();

            Eigen::VectorXcd h_est = receiver_.getEstimatedPathCoefficients();

            double mse_pilot_0 = (h_true_l0 - h_est).squaredNorm();
            double mse_pilot_1 = (h_true_l1 - h_est).squaredNorm();
            totalSquaredError += (mse_pilot_0 + mse_pilot_1) / 2.0;
        }
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.Q_);
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
};

#endif // SIMULATOR_MISC_H
