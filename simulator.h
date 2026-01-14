/*
 * File:   simulator.h
 * Author: Ito
 *
 * Created on 2024/12/20, 18:31
 */

#ifndef SIMULATOR_H
#define SIMULATOR_H
#define _USE_MATH_DEFINES

#include "parameters.h"
#include "channel.h"
#include "transceiver.h"
#include <fstream>    // std::ofstream のために追加
#include <complex>    // std::abs(std::complex) のために追加

class Simulator
{
public:
    Simulator(const SimulationParameters &params)
        : params_(params),
        W_master_(SimulationParameters::generateW(params.K_, params.Q_, params.NUMBER_OF_FFT)),
        channel_(params, W_master_),
        transceiver_(params, W_master_)
    {
    }

    virtual ~Simulator()
    {
    }

    /**
     * 雑音の分散設定
     * @param EbN0dB EbN0 [dB]
     */
    void setNoiseSD(double EbN0dB)
    {
        noiseSD_ = std::sqrt(std::pow(10.0, -0.1 * EbN0dB) / (double)params_.NUMBER_OF_BIT);
    }

    // 試行回数の設定
    //  void setTrialNum(double EbN0dB) {
    //      if(EbN0dB == 0 || EbN0dB == 5 || EbN0dB == 10) {
    //          NUMBER_OF_TRIAL = 10000;
    //      }
    //      else if(EbN0dB == 15 || EbN0dB == 20) {
    //          NUMBER_OF_TRIAL = 100000;
    //      }
    //      else {
    //          NUMBER_OF_TRIAL = 1000000;
    //      }
    //  }

    void setTrialNum(double trialNum)
    {
        NUMBER_OF_TRIAL = trialNum;
    }

    /**
     * ドップラー周波数設定
     * @param f_d ドップラー周波数
     */
    void setDopplerFrequency(double fd_Ts)
    {
        fd_Ts_ = fd_Ts;
    }

    double getAverageIterations() const
    {
        return last_avg_iterations_;
    }

    /**
     * 数値計算実験
     * @return ビット誤り率のシミュレーション値
     */
    double getBER_EM_Simulation()
    {
        int totalErrorCount = 0;
        double total_iterations_sum = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transceiver_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            transceiver_.setY_(channel_.getH(), noiseSD_);
            double trial_avg_iter = transceiver_.equalizeAndDemodulate();
            totalErrorCount += transceiver_.getBitErrorCount();
            total_iterations_sum += trial_avg_iter;
        }
        last_avg_iterations_ = total_iterations_sum / static_cast<double>(NUMBER_OF_TRIAL);
        return (double)totalErrorCount / ((double)NUMBER_OF_TRIAL * (double)params_.NUMBER_OF_BIT * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

    /**
     * チャネル推定のMSEを計算するシミュレーション
     * @return MSEのシミュレーション値
     */
    double getMSE_simulation()
    {
        double totalSquaredError = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transceiver_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            transceiver_.setY_(channel_.getH(), noiseSD_);
            transceiver_.equalizeAndDemodulate(); // この中でH_estが計算・保存される
            totalSquaredError += transceiver_.getMSE();
        }
        // 試行回数、データシンボル数、サブキャリア数で平均化
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

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

    /**
     * 数値計算実験
     * @return パイロットシンボルのみをもちいたビット誤り率のシミュレーション値
     */
    double getBER_Simulation_only_pilot()
    {
        int totalErrorCount = 0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transceiver_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            transceiver_.setY_(channel_.getH(), noiseSD_);
            transceiver_.equalizeByPilotAndDemodulate();
            totalErrorCount += transceiver_.getBitErrorCount();
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
            transceiver_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            transceiver_.setY_(channel_.getH(), noiseSD_);
            transceiver_.equalizeByPilotAndDemodulate(); // この中でH_estが計算・保存される
            totalSquaredError += transceiver_.getMSE();
        }
        // 試行回数、データシンボル数、サブキャリア数で平均化
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

    /**
     * 雑音分散のMSEを計算するシミュレーション
     * @param fixedEbN0dB 固定するEb/N0 [dB]
     * @return 雑音分散MSEのシミュレーション値
     */
    double getNoiseVarianceMSE_simulation(int fixedEbN0dB)
    {
        // setNoiseSD() のロジックを逆算し、真の雑音分散 σ_n^2 を計算
        // σ_n = sqrt(10^(-0.1 * EbN0dB) / NUMBER_OF_BIT)
        // σ_n^2 = 10^(-0.1 * EbN0dB) / NUMBER_OF_BIT
        double true_noise_variance = noiseSD_ * noiseSD_;

        double totalSquaredError = 0.0;
        
        // setNoiseSD(fixedEbN0dB) は main.cpp で事前に呼び出されている必要があります。
        
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transceiver_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            transceiver_.setY_(channel_.getH(), noiseSD_);
            // equalizeAndDemodulate の中で noiseVariance_ (推定値) が更新される
            transceiver_.equalizeAndDemodulate(); 
            
            // 推定値を取得
            double est_noise_variance = transceiver_.getEstimatedNoiseVariance(); //
            
            // MSEの累積: (推定値 - 真値)^2
            double mse_trial = std::pow(est_noise_variance - true_noise_variance, 2.0);
            totalSquaredError += mse_trial;
        }
        
        // 試行回数で平均化
        return totalSquaredError / (double)NUMBER_OF_TRIAL;
    }

    void saveChannelMagnitudeResponseToCSV(std::ofstream& ofs, double fd_Ts)
    {
        // 1. チャネルの生成（1試行のみ）
        channel_.generateFrequencyResponse(fd_Ts);

        // 2. 周波数応答Hを取得
        const auto& H = channel_.getH(); // Hの型は、ユーザーのChannelクラスの実装に依存します
        
        // Hのサイズは params_.K_ (サブキャリア数) x params_.L_ (シンボル数) であると仮定
        int K = params_.K_; // サブキャリア数
        int L = params_.L_; // シンボル数 (全シンボル)

        // CSVヘッダ: k,l,Magnitude
        ofs << "l,k,Magnitude" << std::endl;

        for (int l = 0; l < L; ++l) // 時刻インデックス (0からL-1)
        {
            for (int k = 0; k < K; ++k) // 周波数インデックス (0からK-1)
            {
                double magnitude = std::abs(H(l, k));

                ofs << l << "," << k << "," << magnitude << std::endl;
            }
        }
    }

    double get_h_MSE_Simulation_during_pilot()
    {
        double totalSquaredError = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transceiver_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            transceiver_.setY_(channel_.getH(), noiseSD_);
            transceiver_.est_H_by_initial_h(); // この中でH_estが計算・保存される
            totalSquaredError += transceiver_.getMSE_during_pilot();
        }
        // 試行回数、データシンボル数、サブキャリア数で平均化
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

    double get_H_est_MSE_Simulation_during_pilot()
    {
        double totalSquaredError = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transceiver_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            transceiver_.setY_(channel_.getH(), noiseSD_);
            transceiver_.est_H_by_pilot(); // この中でH_estが計算・保存される
            totalSquaredError += transceiver_.getMSE_during_pilot();
        }
        // 試行回数、データシンボル数、サブキャリア数で平均化
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_);
    }

private:
    SimulationParameters params_;
    Eigen::MatrixXcd W_master_;
    Channel channel_;
    Transceiver transceiver_;

    double noiseSD_;     // 雑音の標準偏差
    int NUMBER_OF_TRIAL; // 試行回数
    double fd_Ts_;         // 正規化ドップラー (f_d * T_s): 単位なし

    double last_avg_iterations_ = 0.0;
};

#endif /* SIMULATOR_H */