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

class Simulator
{
public:
    Simulator(const SimulationParameters &params) : params_(params), channel_(params), transceiver_(params)
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

    /**
     * 数値計算実験
     * @return ビット誤り率のシミュレーション値
     */
    double getBER_EM_Simulation()
    {
        int totalErrorCount = 0;
        for (int tri = 0; tri < NUMBER_OF_TRIAL; tri++)
        {
            transceiver_.setX_();
            channel_.generateFrequencyResponse(fd_Ts_);
            transceiver_.setY_(channel_.getH(), noiseSD_);
            transceiver_.equalizeAndDemodulate(noiseSD_);
            totalErrorCount += transceiver_.getBitErrorCount();
        }
        return (double)totalErrorCount / ((double)NUMBER_OF_TRIAL * (double)params_.NUMBER_OF_BIT * (double)params_.K_ * ((double)params_.L_ - params_.NUMBER_OF_PILOT));
    }

private:
    SimulationParameters params_;
    Channel channel_;
    Transceiver transceiver_;

    double noiseSD_;     // 雑音の標準偏差
    int NUMBER_OF_TRIAL; // 試行回数
    double fd_Ts_;         // 正規化ドップラー (f_d * T_s): 単位なし
};

#endif /* SIMULATOR_H */