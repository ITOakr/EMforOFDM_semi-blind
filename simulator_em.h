#ifndef SIMULATOR_EM_H
#define SIMULATOR_EM_H

#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <atomic>
#include <omp.h>
#include <Eigen/Dense>
#include "simulator_base.h"

class SimulatorEM : public SimulatorBase {
public:
    SimulatorEM(const SimulationParameters &params) : SimulatorBase(params) {}

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
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            double trial_avg_iter = receiver_.equalizeAndDemodulate(transmitter_.getX());
            totalErrorCount += receiver_.getBitErrorCount(transmitter_.getTxData());
            total_iterations_sum += trial_avg_iter;
        }
        last_avg_iterations_ = total_iterations_sum / static_cast<double>(NUMBER_OF_TRIAL);
        return (double)totalErrorCount / ((double)NUMBER_OF_TRIAL * (double)params_.NUMBER_OF_BIT * (double)(params_.enableDataPilots ? (params_.K_ - 4) : params_.K_) * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
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
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.equalizeWithWrapperAIC(transmitter_.getX());
            totalSquaredError += receiver_.getMSE();
        }
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

    /**
     * データ部における判定帰還方式（Decision-Directed）のトラッキング性能（MSE）をシミュレーションする。
     * プリアンブル部で Raghavendra GAIC によりパスモデルと初期伝送路を推定し、
     * データ部では判定帰還方式で毎シンボル伝送路を更新し、純粋なデータキャリアのみでMSEを評価する。
     */
    double get_H_MSE_DecisionDirected_Simulation()
    {
        double totalSquaredError = 0.0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            
            // 1. プリアンブル部: Raghavendra GAIC によるパス選択と初期チャネル推定
            receiver_.est_H_by_initial_h_RaghavendraGAIC(transmitter_.getX());
            
            // 2. データ部: 判定帰還（DD）方式によるトラッキング
            receiver_.est_H_by_DecisionDirected(transmitter_.getX());
            
            // 3. MSE計算 (データ部・データキャリアのみ)
            totalSquaredError += receiver_.getMSE_during_data_only_data_carriers();
        }
        return totalSquaredError / static_cast<double>(NUMBER_OF_TRIAL);
    }


    /**
     * チャネル推定のMSEを計算するシミュレーション (並列化対応版)
     * @return MSEのシミュレーション値
     */
    double getMSE_simulation_parallel()
    {
        double totalSquaredError = 0.0;
        std::atomic<int> completed_trials(0);

        #pragma omp parallel reduction(+:totalSquaredError)
        {
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 

            Channel local_channel(local_params, W_master_);
            Transmitter local_transmitter(local_params);
            Receiver local_receiver(local_params, W_master_);

            #pragma omp for schedule(dynamic)
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transmitter.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_receiver.setY_(local_channel.getH(), local_transmitter.getX(), noiseSD_);
                local_receiver.equalizeWithWrapperAIC(local_transmitter.getX()); 
                totalSquaredError += local_receiver.getMSE();

                int current_count = ++completed_trials;
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
        std::cout << "\rProgress: 100% (" << NUMBER_OF_TRIAL << "/" << NUMBER_OF_TRIAL << ") Done." << std::endl;

        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

    /**
     * 雑音分散のMSEを計算するシミュレーション
     */
    double getNoiseVarianceMSE_simulation(int fixedEbN0dB)
    {
        double true_noise_variance = noiseSD_ * noiseSD_;
        double totalSquaredError = 0.0;
        
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transmitter_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            receiver_.setY_(channel_.getH(), transmitter_.getX(), noiseSD_);
            receiver_.equalizeAndDemodulate(transmitter_.getX()); 
            
            double est_noise_variance = receiver_.getEstimatedNoiseVariance();
            double mse_trial = std::pow(est_noise_variance - true_noise_variance, 2.0);
            totalSquaredError += mse_trial;
        }
        
        return totalSquaredError / (double)NUMBER_OF_TRIAL;
    }

    /**
     * 埋め込み法 (Embedded AIC) によるMSEシミュレーション
     */
    double getMSE_EmbeddedAIC_Simulation()
    {
        double totalSquaredError = 0.0;
        std::atomic<int> completed_trials(0);

        #pragma omp parallel reduction(+:totalSquaredError)
        {
            SimulationParameters local_params = params_;
            local_params.seed += omp_get_thread_num(); 
            Channel local_channel(local_params, W_master_);
            Transmitter local_transmitter(local_params);
            Receiver local_receiver(local_params, W_master_);

            #pragma omp for schedule(dynamic)
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transmitter.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_receiver.setY_(local_channel.getH(), local_transmitter.getX(), noiseSD_);
                local_receiver.equalizeWithEmbeddedAIC(local_transmitter.getX());
                totalSquaredError += local_receiver.getMSE();

                int current_count = ++completed_trials;
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
        std::cout << "\rProgress: 100% (" << NUMBER_OF_TRIAL << "/" << NUMBER_OF_TRIAL << ") Done.   " << std::endl;
        
        return totalSquaredError / ((double)NUMBER_OF_TRIAL * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

    /**
     * Wrapper法 (Wrapper AIC) によるSNR劣化比シミュレーション
     */
    double getSNRDegradation_WrapperAIC_Simulation()
    {
        double totalMetric = 0.0;
        std::atomic<int> completed_trials(0);

        #pragma omp parallel reduction(+:totalMetric)
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
                local_receiver.equalizeWithWrapperAIC(local_transmitter.getX());
                
                totalMetric += local_receiver.getSNRDegradationMetric(noiseSD_);

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
     * パイロットシンボルのみを用いた推定によるSNR劣化比シミュレーション
     */
    double getSNRDegradation_PilotOnly_Simulation()
    {
        double totalMetric = 0.0;
        std::atomic<int> completed_trials(0);

        #pragma omp parallel reduction(+:totalMetric)
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
                local_receiver.equalizeByPilotAndDemodulate(local_transmitter.getX());
                
                totalMetric += local_receiver.getSNRDegradationMetric(noiseSD_);

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
            Transmitter local_transmitter(local_params);
            Receiver local_receiver(local_params, W_master_);

            #pragma omp for
            for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
            {
                local_transmitter.setX_();
                local_channel.generateFrequencyResponse(fd_Ts_);
                local_receiver.setY_(local_channel.getH(), local_transmitter.getX(), noiseSD_);
                local_receiver.equalizeWithPilotAICFixedPath(local_transmitter.getX());
                
                totalSquaredError += local_receiver.getMSE();

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

private:
    double last_avg_iterations_ = 0.0;
};

#endif // SIMULATOR_EM_H
