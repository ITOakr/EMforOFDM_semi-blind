#ifndef SIMULATOR_BASE_H
#define SIMULATOR_BASE_H

#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include "parameters.h"
#include "Channel.h"
#include "transmitter.h"
#include "receiver.h"

class SimulatorBase {
protected:
    SimulationParameters params_;
    Eigen::MatrixXcd W_master_;
    Channel channel_;
    Transmitter transmitter_;
    Receiver receiver_;
    int NUMBER_OF_TRIAL = 1000;
    double noiseSD_ = 0.0;
    double fd_Ts_ = 0.0;

public:
    SimulatorBase(const SimulationParameters &params)
        : params_(params),
          W_master_(SimulationParameters::generateW(params.K_, params.Q_, params.NUMBER_OF_FFT)),
          channel_(params_, W_master_),
          transmitter_(params_),
          receiver_(params_, W_master_) {}

    virtual ~SimulatorBase() {}

    void setNoiseSD(double EbN0dB) {
        noiseSD_ = std::sqrt(std::pow(10.0, -0.1 * EbN0dB) / (double)params_.NUMBER_OF_BIT);
    }

    void setTrialNum(double trialNum) {
        NUMBER_OF_TRIAL = trialNum;
    }

    void setDopplerFrequency(double fd_Ts) {
        fd_Ts_ = fd_Ts;
    }
};

#endif // SIMULATOR_BASE_H
