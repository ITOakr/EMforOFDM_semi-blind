#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include "simulator_em.h"
#include "parameters.h"

// パラメータ
static const int EbN0dBmin = 0;
static const int EbN0dBmax = 30;
static const int EbN0dBstp = 5;

static const int dopplerMin = 0;
static const double dopplerMax = 0.02;
static const double dopplerStep = 0.0002;

std::string fileName;
std::ofstream ofs;
const std::string outputDir = "C:/Users/Akira Ito/code2025/results/Sim_20260204~/";

double ber;
double mse;
int numberOfTrials;

std::string getCurrentTimeString() {
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::tm tm_now;
    
#if defined(_MSC_VER)
    localtime_s(&tm_now, &time_t_now);
#else
    tm_now = *std::localtime(&time_t_now);
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm_now, "%Y%m%d_%H%M%S");
    return oss.str();
}

std::string getModulationSchemeName(int numberOfBits) {
    switch (numberOfBits) {
        case 2: return "QPSK";
        case 4: return "16QAM";
        case 6: return "64QAM";
        case 8: return "256QAM";
        default: return "Unknown";
    }
}

int main()
{
    SimulationParameters params;

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Number of Bit? (QPSK:2, 16QAM:4, 64QAM:6, 256QAM:8)" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> params.NUMBER_OF_BIT;
    params.NUMBER_OF_SYMBOLS = std::pow(2, params.NUMBER_OF_BIT);

    SimulatorEM sim(params);

    int mode_select = 0;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Select EM / Data-Aided Simulation Mode" << std::endl;
    std::cout << "1: BER vs Eb/N0 sweep (fixed Doppler)" << std::endl;
    std::cout << "2: BER vs Doppler sweep (fixed Eb/N0)" << std::endl;
    std::cout << "3: MSE vs Eb/N0 sweep (fixed Doppler)" << std::endl;
    std::cout << "8: MSE vs Doppler sweep (fixed Eb/N0, 0 -> step -> double)" << std::endl;
    std::cout << "11: MSE of estimated noise variance vs Doppler sweep (fixed Eb/N0)" << std::endl;
    std::cout << "14: Parallel MSE vs Doppler sweep (fixed Eb/N0)" << std::endl;
    std::cout << "18: Embedded AIC method MSE vs Doppler (fixed Eb/N0)" << std::endl;
    std::cout << "19: Wrapper AIC method SNR degradation ratio vs Doppler (fixed Eb/N0)" << std::endl;
    std::cout << "20: Pilot-only SNR degradation ratio vs Doppler (fixed Eb/N0)" << std::endl;
    std::cout << "21: Pilot AIC fixed-path MSE vs Doppler (fixed Eb/N0)" << std::endl;
    std::cout << "29: Sweep frame length L and evaluate MSE (fixed Eb/N0 and Doppler)" << std::endl;
    std::cout << "56: MSE of CFR (H) by Decision-Directed tracking (Data carriers only)" << std::endl;
    std::cout << "57: MSE of CFR (H) by Preamble Only (Data carriers only)" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> mode_select;

    std::cout << "Enter number of trials:";
    std::cin >> numberOfTrials;
    sim.setTrialNum(numberOfTrials);

    std::string modulationName = getModulationSchemeName(params.NUMBER_OF_BIT);
    std::string timeStr = getCurrentTimeString();

    if (mode_select == 1)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_BER_vs_EbN0_2path.csv";
        ofs.open(fileName);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            ber = sim.getBER_EM_Simulation();
            double avg_iter = sim.getAverageIterations();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", BER = " << ber << std::endl;
            std::cout << "Avg_Iteration = " << avg_iter << std::endl;
            ofs << EbN0dB << "," << ber << std::endl;
        }
    }
    else if (mode_select == 2)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_BER_vs_Doppler_2path.csv";
        ofs.open(fileName);

        sim.setNoiseSD(fixedEbN0dB);

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; dopplerFrequency += dopplerStep) {
            sim.setDopplerFrequency(dopplerFrequency);

            ber = sim.getBER_EM_Simulation();
            
            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", BER = " << ber << std::endl;
            ofs << dopplerFrequency << "," << ber << std::endl;
        }
    }
    else if (mode_select == 3)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "MSE_2path.csv";
        ofs.open(fileName);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            mse = sim.getMSE_simulation();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
            ofs << EbN0dB << "," << mse << std::endl;
        }
    }
    else if (mode_select == 8)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_MSE_vs_Doppler_2path.csv";
        ofs.open(fileName);

        sim.setNoiseSD(fixedEbN0dB);

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);

            mse = sim.getMSE_simulation();
            
            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", MSE = " << mse << std::endl;
            ofs << dopplerFrequency << "," << mse << std::endl;

            if (dopplerFrequency == 0.0) {
                dopplerFrequency = dopplerStep;
            } else {
                dopplerFrequency *= 2.0;
            }
        }
    }
    else if (mode_select == 11)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_NoiseVar_MSE_vs_Doppler_2path.csv";
        ofs.open(fileName);

        sim.setNoiseSD(fixedEbN0dB);
        double noise_var_mse;

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; dopplerFrequency += dopplerStep) {
            sim.setDopplerFrequency(dopplerFrequency);
            
            noise_var_mse = sim.getNoiseVarianceMSE_simulation(fixedEbN0dB);
            
            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", Noise Var MSE = " << noise_var_mse << std::endl;
            ofs << dopplerFrequency << "," << noise_var_mse << std::endl;
        }
    }
    else if (mode_select == 14)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_MSE_vs_Doppler_2path_parallel.csv";
        ofs.open(fileName);

        sim.setNoiseSD(fixedEbN0dB);

        auto start_total = std::chrono::high_resolution_clock::now();

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);

            mse = sim.getMSE_simulation_parallel();

            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", MSE = " << mse << std::endl;
            ofs << dopplerFrequency << "," << mse << std::endl;

            if (dopplerFrequency == 0.0) {
                dopplerFrequency = dopplerStep;
            } else {
                dopplerFrequency *= 2.0;
            }
        }
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        
        std::cout << "========================================" << std::endl;
        std::cout << "Total Simulation Time: " << elapsed_total.count() << " seconds." << std::endl;
        std::cout << "========================================" << std::endl;
    }
    else if (mode_select == 18)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "_EmbeddedAIC_MSE_vs_Doppler_EbN0_" + std::to_string(fixedEbN0dB) + ".csv";
        ofs.open(fileName);

        sim.setNoiseSD((double)fixedEbN0dB);

        auto start_total = std::chrono::high_resolution_clock::now();

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);

            mse = sim.getMSE_EmbeddedAIC_Simulation();

            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", MSE = " << mse << std::endl;
            ofs << dopplerFrequency << "," << mse << std::endl;

            if (dopplerFrequency == 0.0) dopplerFrequency = dopplerStep;
            else dopplerFrequency *= 2.0;
        }
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        std::cout << "Total Time: " << elapsed_total.count() << "s" << std::endl;
    }
    else if (mode_select == 19)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "_WrapperAIC_SNRDegradation_vs_Doppler_EbN0_" + std::to_string(fixedEbN0dB) + ".csv";
        ofs.open(fileName);

        sim.setNoiseSD((double)fixedEbN0dB);

        auto start_total = std::chrono::high_resolution_clock::now();

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);
            std::cout << "Target: f_dT_s = " << dopplerFrequency << " processing..." << std::endl;

            double degradation = sim.getSNRDegradation_WrapperAIC_Simulation();

            std::cout << " Result: SNR degradation ratio = " << degradation << std::endl;
            ofs << dopplerFrequency << "," << degradation << std::endl;

            if (dopplerFrequency == 0.0) dopplerFrequency = dopplerStep;
            else dopplerFrequency *= 2.0;
        }
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        std::cout << "Total Time: " << elapsed_total.count() << "s" << std::endl;
    }
    else if (mode_select == 20)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "_PilotOnly_SNRDegradation_vs_Doppler_EbN0_" + std::to_string(fixedEbN0dB) + ".csv";
        ofs.open(fileName);

        sim.setNoiseSD((double)fixedEbN0dB);

        auto start_total = std::chrono::high_resolution_clock::now();

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);
            std::cout << "Target: f_dT_s = " << dopplerFrequency << " processing..." << std::endl;

            double degradation = sim.getSNRDegradation_PilotOnly_Simulation();

            std::cout << " Result: SNR degradation ratio = " << degradation << std::endl;
            ofs << dopplerFrequency << "," << degradation << std::endl;

            if (dopplerFrequency == 0.0) dopplerFrequency = dopplerStep;
            else dopplerFrequency *= 2.0;
        }
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        std::cout << "Total Time: " << elapsed_total.count() << "s" << std::endl;
    }
    else if (mode_select == 21)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "_PilotAICFixed_MSE.csv";
        ofs.open(fileName);

        sim.setNoiseSD((double)fixedEbN0dB);
        auto start_total = std::chrono::high_resolution_clock::now();

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);
            std::cout << "Target: f_dT_s = " << dopplerFrequency << " processing..." << std::endl;

            double mse_val = sim.getMSE_PilotAICFixedPath_Simulation();

            std::cout << " Result: MSE = " << mse_val << std::endl;
            ofs << dopplerFrequency << "," << mse_val << std::endl;

            if (dopplerFrequency == 0.0) dopplerFrequency = dopplerStep;
            else dopplerFrequency *= 2.0;
        }
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        std::cout << "Total Time: " << elapsed_total.count() << "s" << std::endl;
    }
    else if (mode_select == 29)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:";
        std::cin >> fixedEbN0dB;

        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        int L_start, L_end, L_step;
        std::cout << "Enter L start (min 2): "; std::cin >> L_start;
        std::cout << "Enter L end: "; std::cin >> L_end;
        std::cout << "Enter L step: "; std::cin >> L_step;

        fileName = outputDir + timeStr + "_" + modulationName + "_EbN0_" + std::to_string(fixedEbN0dB) + "_fdTs_" + std::to_string(dopplerFrequency) + "_MSE_vs_L.csv";
        ofs.open(fileName);
        ofs << "L,MSE" << std::endl;

        std::cout << "Starting L sweep simulation..." << std::endl;
        auto start_total = std::chrono::high_resolution_clock::now();

        for (int L = L_start; L <= L_end; L += L_step) {
            params.L_ = L;
            
            SimulatorEM sim_L(params);
            sim_L.setTrialNum(numberOfTrials);
            sim_L.setDopplerFrequency(dopplerFrequency);
            sim_L.setNoiseSD(fixedEbN0dB);

            double mse_val = sim_L.getMSE_simulation_parallel();

            std::cout << "L = " << L << ", MSE = " << mse_val << std::endl;
            ofs << L << "," << mse_val << std::endl;
        }

        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_total = end_total - start_total;
        std::cout << "Total Simulation Time: " << elapsed_total.count() << " seconds." << std::endl;
        std::cout << "Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 56)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_DecisionDirected_DataCarrierMSE_Mode56.csv";
        ofs.open(fileName);
        std::cout << "Starting Decision-Directed Tracking MSE Simulation (Mode 56)..." << std::endl;
        std::cout << "EbN0dB, MSE" << std::endl;
        ofs << "EbN0dB,MSE" << std::endl;

        sim.setDopplerFrequency(dopplerFrequency);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setNoiseSD(EbN0dB);
            double mse_val = sim.get_H_MSE_DecisionDirected_Simulation();
            std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse_val << std::endl;
            ofs << EbN0dB << "," << mse_val << std::endl;
        }
        std::cout << "Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 57)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_PreambleOnly_DataCarrierMSE_Mode57.csv";
        ofs.open(fileName);
        std::cout << "Starting Preamble Only Tracking MSE Simulation (Mode 57)..." << std::endl;
        std::cout << "EbN0dB, MSE" << std::endl;
        ofs << "EbN0dB,MSE" << std::endl;

        sim.setDopplerFrequency(dopplerFrequency);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setNoiseSD(EbN0dB);
            double mse_val = sim.get_H_MSE_PreambleOnly_Simulation();
            std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse_val << std::endl;
            ofs << EbN0dB << "," << mse_val << std::endl;
        }
        std::cout << "Results saved to: " << fileName << std::endl;
    }
    else
    {
        std::cout << "Invalid mode selected." << std::endl;
        if (ofs.is_open()) ofs.close();
        return 1;
    }

    if (ofs.is_open()) ofs.close();
    return 0;
}
