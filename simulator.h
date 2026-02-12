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
#include <omp.h>      // OpenMP のために追加
#include <atomic> // 追加: 進捗カウント用

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
            transceiver_.equalizeWithWrapperAIC(); // この中でH_estが計算・保存される(equalizeAndDemodulate())
            totalSquaredError += transceiver_.getMSE();
        }
        // 試行回数、データシンボル数、サブキャリア数で平均化
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

    /**
     * チャネル推定のMSEを計算するシミュレーション (並列化対応版)
     * @return MSEのシミュレーション値
     */
    double getMSE_simulation_parallel()
    {
        double totalSquaredError = 0.0;
        
        // ★進捗カウント用の変数を定義（並列ブロックの外で作る）
        // std::atomic を使うことで、複数のスレッドから同時に +1 されても正しくカウントされます
        std::atomic<int> completed_trials(0);

        // OpenMPによる並列化
        // reduction(+:totalSquaredError) で各スレッドの結果を合計する
        #pragma omp parallel reduction(+:totalSquaredError)
        {
            // --- スレッドローカルな変数の準備 ---
            
            // 1. パラメータをコピーし、乱数シードをスレッドID分ずらす
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 

            // 2. このスレッド専用の Channel と Transceiver を作成
            // (メンバ変数の channel_, transceiver_ は使わず、ここ専用のものを使う)
            // W_master_ は読み取り専用なので共有のものを渡してOK
            Channel local_channel(local_params, W_master_);
            Transceiver local_transceiver(local_params, W_master_);

            // --- ループ処理 ---
            #pragma omp for schedule(dynamic)
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transceiver.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                
                // noiseSD_ は共有変数だが読み取り専用なのでそのままでOK
                local_transceiver.setY_(local_channel.getH(), noiseSD_);
                
                // Wrapper法の実行
                local_transceiver.equalizeWithWrapperAIC(); 
                
                // MSEの蓄積
                totalSquaredError += local_transceiver.getMSE();

                // 1. 完了数を +1 する (スレッドセーフ)
                int current_count = ++completed_trials;

                // 2. 一定間隔（例: 10%ごと）でのみ表示処理を行う
                // 頻繁に cout すると逆に遅くなるため、間引きが必要です
                if (NUMBER_OF_TRIAL >= 10 && (current_count % (NUMBER_OF_TRIAL / 10) == 0)) 
                {
                    // 3. 表示中は他のスレッドを待たせる（文字化け防止）
                    #pragma omp critical
                    {
                        // \r で行頭に戻って上書き表示する
                        double progress = (double)current_count / NUMBER_OF_TRIAL * 100.0;
                        std::cout << "\rProgress: " << (int)progress << "% (" << current_count << "/" << NUMBER_OF_TRIAL << ")" << std::flush;
                    }
                }
            }
        } // 並列領域終了
        std::cout << "\rProgress: 100% (" << NUMBER_OF_TRIAL << "/" << NUMBER_OF_TRIAL << ") Done." << std::endl;

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

    /**
     * AICによるモデル選択の正答率を計算するシミュレーション
     * @return 正答率 (0.0 ~ 1.0)
     */
    double getAICAccuracy_pilot()
    {
        int successCount = 0;
        
        // 並列化する場合 (高速化のため推奨)
        // std::atomic<int> completed_trials(0); // 進捗表示用
        #pragma omp parallel reduction(+:successCount)
        {
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 
            Channel local_channel(local_params, W_master_);
            Transceiver local_transceiver(local_params, W_master_);

            #pragma omp for
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transceiver.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_transceiver.setY_(local_channel.getH(), noiseSD_);

                // パイロットのみでAIC推定を実行 (set_initial_params_by_pilotが内部で呼ばれる)
                local_transceiver.est_H_by_initial_h(); 

                // 正解判定 (pathMaskと比較)
                if (local_transceiver.checkAICAccuracy(local_params.pathMask)) {
                    successCount++;
                }
            }
        }

        return (double)successCount / (double)NUMBER_OF_TRIAL;
    }

    /**
     * AICの評価指標（F値と正答率）を計算する
     * @return pair<F-measure, Accuracy>
     */
    std::pair<double, double> getAIC_Metrics_pilot()
    {
        long long total_TP = 0; // 正解: あるものを「ある」と言えた
        long long total_FP = 0; // 誤検知: ないものを「ある」と言ってしまった
        long long total_FN = 0; // 見逃し: あるものを「ない」と言ってしまった
        long long total_TN = 0; // 正常: ないものを「ない」と言えた

        #pragma omp parallel reduction(+:total_TP, total_FP, total_FN, total_TN)
        {
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 
            Channel local_channel(local_params, W_master_);
            Transceiver local_transceiver(local_params, W_master_);

            #pragma omp for
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transceiver.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_transceiver.setY_(local_channel.getH(), noiseSD_);

                local_transceiver.est_H_by_initial_h(); 
                std::vector<int> estMask = local_transceiver.getEstimatedPathMask();
                
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

        // --- F値の計算 ---
        long long denominator_prec = total_TP + total_FP;
        double precision = (denominator_prec > 0) ? (double)total_TP / denominator_prec : 0.0;

        long long denominator_rec = total_TP + total_FN;
        double recall = (denominator_rec > 0) ? (double)total_TP / denominator_rec : 0.0;

        double numerator_f = 2.0 * precision * recall;
        double denominator_f = precision + recall;
        double f_measure = (denominator_f > 0) ? numerator_f / denominator_f : 0.0;

        // --- 正答率の計算 ---
        long long total = total_TP + total_FP + total_FN + total_TN;
        double accuracy = (total > 0) ? (double)(total_TP + total_TN) / total : 0.0;

        return {f_measure, accuracy};
    }

    /**
     * 埋め込み法 (Embedded AIC) によるMSEシミュレーション (進捗表示あり)
     * @return MSE
     */
    double getMSE_EmbeddedAIC_Simulation()
    {
        double totalSquaredError = 0.0;
        
        // 進捗カウント用 (アトミック変数)
        std::atomic<int> completed_trials(0);

        #pragma omp parallel reduction(+:totalSquaredError)
        {
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 
            Channel local_channel(local_params, W_master_);
            Transceiver local_transceiver(local_params, W_master_);

            #pragma omp for schedule(dynamic)
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transceiver.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_transceiver.setY_(local_channel.getH(), noiseSD_);
                
                // 埋め込み法を実行
                local_transceiver.equalizeWithEmbeddedAIC();
                
                totalSquaredError += local_transceiver.getMSE();

                // --- 進捗表示ロジック ---
                int current_count = ++completed_trials;
                // 10%刻み、または100回に1回など適度な頻度で表示
                if (NUMBER_OF_TRIAL >= 10 && (current_count % (NUMBER_OF_TRIAL / 10) == 0)) 
                {
                    #pragma omp critical
                    {
                        double progress = (double)current_count / NUMBER_OF_TRIAL * 100.0;
                        std::cout << "\rProgress: " << (int)progress << "% (" << current_count << "/" << NUMBER_OF_TRIAL << ")" << std::flush;
                    }
                }
            }
        }
        // ループ終了後に改行
        std::cout << "\rProgress: 100% (" << NUMBER_OF_TRIAL << "/" << NUMBER_OF_TRIAL << ") Done.   " << std::endl;
        
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

    /**
     * Wrapper法 (Wrapper AIC) によるSNR劣化比シミュレーション (並列化・進捗表示)
     * @return 平均SNR劣化比
     */
    double getSNRDegradation_WrapperAIC_Simulation()
    {
        double totalMetric = 0.0;
        
        // 進捗カウント用の変数
        std::atomic<int> completed_trials(0);

        #pragma omp parallel reduction(+:totalMetric)
        {
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 
            Channel local_channel(local_params, W_master_);
            Transceiver local_transceiver(local_params, W_master_);

            #pragma omp for
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transceiver.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_transceiver.setY_(local_channel.getH(), noiseSD_);
                
                // Wrapper法を実行
                local_transceiver.equalizeWithWrapperAIC();
                
                // 評価指標の計算 (noiseSDが必要)
                totalMetric += local_transceiver.getSNRDegradationMetric(noiseSD_);

                // --- 進捗表示ロジック ---
                int current_count = ++completed_trials;
                if (NUMBER_OF_TRIAL >= 10 && (current_count % (NUMBER_OF_TRIAL / 10) == 0)) 
                {
                    #pragma omp critical
                    {
                        double progress = (double)current_count / NUMBER_OF_TRIAL * 100.0;
                        std::cout << "\rProgress: " << (int)progress << "% (" << current_count << "/" << NUMBER_OF_TRIAL << ")   " << std::flush;
                    }
                }
            }
        }
        
        std::cout << "\rProgress: 100% (" << NUMBER_OF_TRIAL << "/" << NUMBER_OF_TRIAL << ") Done.   " << std::endl;
        
        return totalMetric / (double)NUMBER_OF_TRIAL;
    }

    /**
     * パイロットシンボルのみを用いた推定によるSNR劣化比シミュレーション (並列化・進捗表示)
     * @return 平均SNR劣化比
     */
    double getSNRDegradation_PilotOnly_Simulation()
    {
        double totalMetric = 0.0;
        
        // 進捗カウント用の変数 (スレッドセーフ)
        std::atomic<int> completed_trials(0);

        #pragma omp parallel reduction(+:totalMetric)
        {
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 
            Channel local_channel(local_params, W_master_);
            Transceiver local_transceiver(local_params, W_master_);

            #pragma omp for
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transceiver.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_transceiver.setY_(local_channel.getH(), noiseSD_);
                
                // パイロットシンボルのみで等化（AIC推定含む）
                // これにより H_est_ がパイロット時点の推定値で埋められます
                local_transceiver.equalizeByPilotAndDemodulate();
                
                // 評価指標の計算
                totalMetric += local_transceiver.getSNRDegradationMetric(noiseSD_);

                // --- 進捗表示ロジック ---
                int current_count = ++completed_trials;
                // 10%刻みで表示
                if (NUMBER_OF_TRIAL >= 10 && (current_count % (NUMBER_OF_TRIAL / 10) == 0)) 
                {
                    #pragma omp critical
                    {
                        double progress = (double)current_count / NUMBER_OF_TRIAL * 100.0;
                        std::cout << "\rProgress: " << (int)progress << "% (" << current_count << "/" << NUMBER_OF_TRIAL << ")   " << std::flush;
                    }
                }
            }
        }
        
        std::cout << "\rProgress: 100% (" << NUMBER_OF_TRIAL << "/" << NUMBER_OF_TRIAL << ") Done.   " << std::endl;
        
        return totalMetric / (double)NUMBER_OF_TRIAL;
    }

    /**
     * パイロットAIC固定パス法によるMSEシミュレーション
     */
    double getMSE_PilotAICFixedPath_Simulation()
    {
        double totalSquaredError = 0.0;
        std::atomic<int> completed_trials(0);

        #pragma omp parallel reduction(+:totalSquaredError)
        {
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 
            Channel local_channel(local_params, W_master_);
            Transceiver local_transceiver(local_params, W_master_);

            #pragma omp for
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transceiver.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_transceiver.setY_(local_channel.getH(), noiseSD_);
                
                // 新しい手法を実行
                local_transceiver.equalizeWithPilotAICFixedPath();
                
                totalSquaredError += local_transceiver.getMSE();

                // 進捗表示
                int current_count = ++completed_trials;
                if (NUMBER_OF_TRIAL >= 10 && (current_count % (NUMBER_OF_TRIAL / 10) == 0)) 
                {
                    #pragma omp critical
                    {
                        double progress = (double)current_count / NUMBER_OF_TRIAL * 100.0;
                        std::cout << "\rProgress: " << (int)progress << "% (" << current_count << "/" << NUMBER_OF_TRIAL << ")   " << std::flush;
                    }
                }
            }
        }
        std::cout << "\rProgress: 100% (" << NUMBER_OF_TRIAL << "/" << NUMBER_OF_TRIAL << ") Done.   " << std::endl;

        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

    /**
     * [Mode 23] 送信信号波形の出力実行関数
     */
    void runExportTxWaveform(int target_k, const std::string& filename)
    {
        // 1回だけ試行環境を作る
        SimulationParameters local_params = params_;
        Transceiver local_transceiver(local_params, W_master_);

        // 信号生成
        local_transceiver.setX_();
        
        // 出力
        local_transceiver.exportTxSymbolTrace(target_k, filename);
    }

    /**
     * [Mode 24] フェージングを受けた信号波形の出力実行関数
     */
    void runExportFadedWaveform(int target_k, const std::string& filename)
    {
        // 1回だけ試行環境を作る
        SimulationParameters local_params = params_;
        Channel local_channel(local_params, W_master_);
        Transceiver local_transceiver(local_params, W_master_);

        // 信号生成
        local_transceiver.setX_();
        // チャネル生成 (ここでドップラー周波数 fd_Ts_ が使われます)
        local_channel.generateFrequencyResponse(fd_Ts_);
        
        // 出力 (生成されたHを渡す)
        local_transceiver.exportFadedSymbolTrace(target_k, filename, local_channel.getH());
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