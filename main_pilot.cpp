#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include "simulator_pilot.h"
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

    SimulatorPilot sim(params);

    int mode_select = 0;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Select Pilot-Only Simulation Mode" << std::endl;
    std::cout << "5: BER by Basic Pilot LS Estimation vs Eb/N0" << std::endl;
    std::cout << "6: MSE of CFR (H) by Basic Pilot LS Estimation vs Eb/N0" << std::endl;
    std::cout << "7: BER by Basic Pilot LS Estimation vs Doppler" << std::endl;
    std::cout << "9: MSE of CFR (H) by Basic Pilot LS Estimation vs Doppler" << std::endl;
    std::cout << "12: MSE of CIR (h) by Full Q-path LS Estimation vs Eb/N0" << std::endl;
    std::cout << "13: MSE of CFR (H) by Full Q-path LS Estimation vs Eb/N0" << std::endl;
    std::cout << "15: Accuracy of Standard AIC Path Selection vs Doppler" << std::endl;
    std::cout << "16: Accuracy of Path Selection (Metric 2) vs Doppler" << std::endl;
    std::cout << "17: F-measure of Standard AIC Path Selection vs Doppler" << std::endl;
    std::cout << "37: MSE of CFR (H) by Fixed 16-path LS Estimation vs Eb/N0" << std::endl;
    std::cout << "38: MSE of CFR (H) by Raghavendra AIC vs Eb/N0" << std::endl;
    std::cout << "39: Export Average CIR Power by Raghavendra AIC" << std::endl;
    std::cout << "40: MSE of CIR (h) by Full LS Estimation (Random Paths) vs Eb/N0" << std::endl;
    std::cout << "41: Accuracy of Exhaustive AIC (8 paths) Path Selection vs Eb/N0" << std::endl;
    std::cout << "42: MSE of CFR (H) by Raghavendra AIC (Random Paths) vs Eb/N0" << std::endl;
    std::cout << "43: MSE of CFR (H) by Exhaustive AIC (8 paths, Fixed Mask) vs Eb/N0" << std::endl;
    std::cout << "44: MSE of CFR (H) by Exhaustive Raghavendra AIC (8 paths, Fixed Mask) vs Eb/N0" << std::endl;
    std::cout << "45: MSE of CFR (H) by Raghavendra AIC Update (Fixed Mask) vs Eb/N0" << std::endl;
    std::cout << "46: MSE of CFR (H) by Known Mask LS Estimation (Random Paths) vs Eb/N0" << std::endl;
    std::cout << "47: Comparison: Standard AIC vs Raghavendra GAIC (Single Trial)" << std::endl;
    std::cout << "48: Comparison: Standard AIC vs Raghavendra GAIC (Average over Trials)" << std::endl;
    std::cout << "49: MSE of CFR (H) by Power-sort Raghavendra GAIC vs Eb/N0" << std::endl;
    std::cout << "50: Comparison: Standard AIC vs Raghavendra GAIC (Fixed Mask, Average over Trials)" << std::endl;
    std::cout << "51: MSE of CFR (H) by Simplified Power-sort Raghavendra GAIC vs Eb/N0" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> mode_select;

    std::cout << "Enter number of trials:";
    std::cin >> numberOfTrials;
    sim.setTrialNum(numberOfTrials);

    std::string modulationName = getModulationSchemeName(params.NUMBER_OF_BIT);
    std::string timeStr = getCurrentTimeString();

    if (mode_select == 5)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_BER_vs_EbN0_2path_pilot.csv";
        ofs.open(fileName);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            ber = sim.getBER_Simulation_only_pilot();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", BER = " << ber << std::endl;
            ofs << EbN0dB << "," << ber << std::endl;
        }
    }
    else if (mode_select == 6)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "MSE_pilot.csv";
        ofs.open(fileName);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            mse = sim.getMSE_simulation_only_pilot();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
            ofs << EbN0dB << "," << mse << std::endl;
        }
    }
    else if (mode_select == 7)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_BER_vs_Doppler_2path_only_pilot.csv";
        ofs.open(fileName);

        sim.setNoiseSD(fixedEbN0dB);

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; dopplerFrequency += dopplerStep) {
            sim.setDopplerFrequency(dopplerFrequency);

            ber = sim.getBER_Simulation_only_pilot();
            
            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", BER = " << ber << std::endl;
            ofs << dopplerFrequency << "," << ber << std::endl;
        }
    }
    else if (mode_select == 9)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:" << std::endl;
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_MSE_vs_Doppler_2path_only_pilot.csv";
        ofs.open(fileName);

        sim.setNoiseSD(fixedEbN0dB);

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);

            mse = sim.getMSE_simulation_only_pilot();
            
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
    else if (mode_select == 12)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_MSE_vs_EbN0_pilot_h_est_MODE12.csv";
        ofs.open(fileName);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            mse = sim.get_h_MSE_Simulation_during_pilot();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
            ofs << EbN0dB << "," << mse << std::endl;
        }
    }
    else if (mode_select == 13)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_MSE_vs_EbN0_pilot_H_est_MODE13.csv";
        ofs.open(fileName);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            mse = sim.get_H_est_MSE_Simulation_during_pilot();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
            ofs << EbN0dB << "," << mse << std::endl;
        }
    }
    else if (mode_select == 15)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_AIC_Model_Accuracy_vs_Doppler.csv";
        ofs.open(fileName);

        sim.setNoiseSD(fixedEbN0dB);
        double accuracy;

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);

            accuracy = sim.getAICAccuracy_pilot();
            
            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", Accuracy = " << accuracy << std::endl;
            ofs << dopplerFrequency << "," << accuracy << std::endl;

            if (dopplerFrequency == 0.0) {
                dopplerFrequency = dopplerStep;
            } else {
                dopplerFrequency *= 2.0;
            }
        }
    }
    else if (mode_select == 16)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_AIC_Path_Accuracy.csv";
        ofs.open(fileName);
        
        ofs << "f_dT_s,Accuracy" << std::endl;

        sim.setNoiseSD(fixedEbN0dB);

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);

            std::pair<double, double> result = sim.getAIC_Metrics_pilot();
            double accuracy = result.second;
            
            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", Accuracy = " << accuracy << std::endl;
            ofs << dopplerFrequency << "," << accuracy << std::endl;

            if (dopplerFrequency == 0.0) dopplerFrequency = dopplerStep;
            else dopplerFrequency *= 2.0;
        }
    }
    else if (mode_select == 17)
    {
        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "EbN0_" + std::to_string(fixedEbN0dB) + "_AIC_F_Measure.csv";
        ofs.open(fileName);
        
        ofs << "f_dT_s,F_Measure" << std::endl;

        sim.setNoiseSD(fixedEbN0dB);

        for (double dopplerFrequency = dopplerMin; dopplerFrequency <= dopplerMax; ) {
            sim.setDopplerFrequency(dopplerFrequency);

            std::pair<double, double> result = sim.getAIC_Metrics_pilot();
            double f_measure = result.first;
            
            std::cout << "-----------" << std::endl;
            std::cout << "f_dT_s = " << dopplerFrequency << ", F-Measure = " << f_measure << std::endl;
            ofs << dopplerFrequency << "," << f_measure << std::endl;

            if (dopplerFrequency == 0.0) dopplerFrequency = dopplerStep;
            else dopplerFrequency *= 2.0;
        }
    }
    else if (mode_select == 37)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_16Paths_HMSE_vs_EbN0.csv";
        ofs.open(fileName);
        ofs << "EbN0dB,H_MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double H_mse = sim.get_H_MSE_with_16paths_during_pilot();

            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", H_MSE = " << H_mse << std::endl;
            ofs << EbN0dB << "," << H_mse << std::endl;
        }
    }
    else if (mode_select == 38)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_RaghavendraAIC_MSE_vs_EbN0.csv";
        ofs.open(fileName);
        ofs << "EbN0dB,H_MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double H_mse = sim.get_H_MSE_with_RaghavendraAIC_during_pilot();

            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", H_MSE = " << H_mse << std::endl;
            ofs << EbN0dB << "," << H_mse << std::endl;
        }
    }
    else if (mode_select == 39)
    {
        double inputDoppler;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> inputDoppler;
        sim.setDopplerFrequency(inputDoppler);

        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:";
        std::cin >> fixedEbN0dB;
        sim.setNoiseSD(fixedEbN0dB);

        fileName = outputDir + timeStr + "_" + modulationName + "_EbN0_" + std::to_string(fixedEbN0dB) + "_16Paths_AverageImpulseResp_q.csv";
        
        std::cout << "Exporting Average Estimated Impulse Response Power (16 paths, l=0)..." << std::endl;
        sim.saveAverageEstimatedImpulseResponseByQ_16paths(fileName);
        
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Saved to: " << fileName << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
    }
    else if (mode_select == 40)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_RandomPath_MSE_MODE40.csv";
        ofs.open(fileName);

        std::cout << "Starting Random Path MSE Simulation (Mode 40)..." << std::endl;
        std::cout << "Eb/N0 [dB], MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double mse_val = sim.getMSE_RandomPath_Mode12_Simulation();

            std::cout << EbN0dB << ", " << mse_val << std::endl;
            ofs << EbN0dB << "," << mse_val << std::endl;
        }
        std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 41)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_ExhaustiveAIC8_Accuracy_EbN0.csv";
        ofs.open(fileName);

        std::cout << "Starting Exhaustive AIC (8 paths) Accuracy Simulation (Mode 41)..." << std::endl;
        std::cout << "Eb/N0 [dB], Accuracy" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double accuracy = sim.getExhaustiveAICAccuracy_8paths_Simulation();

            std::cout << EbN0dB << ", " << accuracy << std::endl;
            ofs << EbN0dB << "," << accuracy << std::endl;
        }
        std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 42)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_RandomPath_MSE_MODE42_Raghavendra.csv";
        ofs.open(fileName);

        std::cout << "Starting Random Path MSE Simulation with Raghavendra AIC (Mode 42)..." << std::endl;
        std::cout << "Eb/N0 [dB], MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double mse_val = sim.getMSE_RandomPath_RaghavendraAIC_Simulation();

            std::cout << EbN0dB << ", " << mse_val << std::endl;
            ofs << EbN0dB << "," << mse_val << std::endl;
        }
        std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 43)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_ExhaustiveAIC8_FixedMask_MSE_EbN0.csv";
        ofs.open(fileName);

        std::cout << "Starting Exhaustive AIC (8 paths) Fixed Mask MSE Simulation (Mode 43)..." << std::endl;
        std::cout << "Eb/N0 [dB], MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double mse_val = sim.getMSE_ExhaustiveAIC_8paths_fixedMask_Simulation();

            std::cout << EbN0dB << ", " << mse_val << std::endl;
            ofs << EbN0dB << "," << mse_val << std::endl;
        }
        std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 44)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_ExhaustiveRaghavendraAIC8_FixedMask_MSE_EbN0.csv";
        ofs.open(fileName);

        std::cout << "Starting Exhaustive Raghavendra AIC (8 paths) Fixed Mask MSE Simulation (Mode 44)..." << std::endl;
        std::cout << "Eb/N0 [dB], MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double mse_val = sim.getMSE_ExhaustiveRaghavendraAIC_8paths_fixedMask_Simulation();

            std::cout << EbN0dB << ", " << mse_val << std::endl;
            ofs << EbN0dB << "," << mse_val << std::endl;
        }
        std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 45)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_RandomPath_MSE_MODE45_Raghavendra.csv";
        ofs.open(fileName);

        std::cout << "Starting Random Path MSE Simulation with Raghavendra AIC (Mode 45)..." << std::endl;
        std::cout << "Eb/N0 [dB], MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double mse_val = sim.getMSE_RandomPath_RaghavendraAIC_Simulation2();

            std::cout << EbN0dB << ", " << mse_val << std::endl;
            ofs << EbN0dB << "," << mse_val << std::endl;
        }
        std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 46)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_RandomPath_MSE_MODE46_KnownMask.csv";
        ofs.open(fileName);

        std::cout << "Starting Random Path MSE Simulation with Known Mask (Mode 46)..." << std::endl;
        std::cout << "Eb/N0 [dB], MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double mse_val = sim.getMSE_RandomPath_KnownMask_Simulation();

            std::cout << EbN0dB << ", " << mse_val << std::endl;
            ofs << EbN0dB << "," << mse_val << std::endl;
        }
        std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 47)
    {
        double dopplerFrequency, fixedEbN0dB;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> dopplerFrequency;
        std::cout << "Enter fixed Eb/N0 [dB] for this trial:";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_EbN0=" + std::to_string(fixedEbN0dB) + "_AIC_vs_GAIC_SingleTrial_MODE47.csv";
        ofs.open(fileName);

        std::cout << "Running AIC vs Raghavendra GAIC Single Trial Simulation (Mode 47)..." << std::endl;
        std::cout << "Q, AIC, Raghavendra_GAIC" << std::endl;
        ofs << "Q,AIC,Raghavendra_GAIC" << std::endl;

        auto [aic, gaic] = sim.getAICvsQ_SingleTrial_Simulation(dopplerFrequency, fixedEbN0dB);

        for (size_t q = 0; q < aic.size(); ++q) {
            std::cout << (q + 1) << ", " << aic[q] << ", " << gaic[q] << std::endl;
            ofs << (q + 1) << "," << aic[q] << "," << gaic[q] << std::endl;
        }
        std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 48)
    {
        double dopplerFrequency, fixedEbN0dB;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> dopplerFrequency;
        std::cout << "Enter fixed Eb/N0 [dB] for average simulation:";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_EbN0=" + std::to_string(fixedEbN0dB) + "_AIC_vs_GAIC_Average_MODE48.csv";
        ofs.open(fileName);

        std::cout << "Running AIC vs Raghavendra GAIC Average Simulation (" << numberOfTrials << " trials) (Mode 48)..." << std::endl;
        std::cout << "Q, Average_AIC, Average_Raghavendra_GAIC" << std::endl;
        ofs << "Q,Average_AIC,Average_Raghavendra_GAIC" << std::endl;

        auto [avg_aic, avg_gaic] = sim.getAICvsQ_Average_Simulation(dopplerFrequency, fixedEbN0dB);

        for (size_t q = 0; q < avg_aic.size(); ++q) {
            std::cout << (q + 1) << ", " << avg_aic[q] << ", " << avg_gaic[q] << std::endl;
            ofs << (q + 1) << "," << avg_aic[q] << "," << avg_gaic[q] << std::endl;
        }
        std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 49)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_MSE_vs_EbN0_pilot_h_est_MODE49_RaghavendraGAIC.csv";
        ofs.open(fileName);

        auto start = std::chrono::high_resolution_clock::now();

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            mse = sim.get_H_MSE_by_pilot_power_sort_RaghavendraGAIC();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
            ofs << EbN0dB << "," << mse << std::endl;
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Mode 49 Total Execution Time: " << elapsed.count() << " seconds" << std::endl;
    }
    else if (mode_select == 50)
    {
        double dopplerFrequency, fixedEbN0dB;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> dopplerFrequency;
        std::cout << "Enter fixed Eb/N0 [dB] for average simulation:";
        std::cin >> fixedEbN0dB;

        fileName = outputDir + timeStr + "_" + modulationName + "_f_dT_s=" + std::to_string(dopplerFrequency) + "_EbN0=" + std::to_string(fixedEbN0dB) + "_AIC_vs_GAIC_Average_FixedMask_MODE50.csv";
        ofs.open(fileName);

        std::cout << "Running AIC vs Raghavendra GAIC Average Simulation (" << numberOfTrials << " trials) (Mode 50)..." << std::endl;
        std::cout << "Q, Average_AIC, Average_Raghavendra_GAIC" << std::endl;
        ofs << "Q,Average_AIC,Average_Raghavendra_GAIC" << std::endl;

        auto [avg_aic, avg_gaic] = sim.getAICvsQ_Average_FixedMask_Simulation(dopplerFrequency, fixedEbN0dB);

        for (size_t q = 0; q < avg_aic.size(); ++q) {
            std::cout << (q + 1) << ", " << avg_aic[q] << ", " << avg_gaic[q] << std::endl;
            ofs << (q + 1) << "," << avg_aic[q] << "," << avg_gaic[q] << std::endl;
        }
        std::cout << "Simulation Completed. Results saved to: " << fileName << std::endl;
    }
    else if (mode_select == 51)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_MSE_vs_EbN0_pilot_h_est_MODE51_RaghavendraGAIC_simplified.csv";
        ofs.open(fileName);

        auto start = std::chrono::high_resolution_clock::now();

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            mse = sim.get_H_MSE_by_pilot_power_sort_RaghavendraGAIC_simplified();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", MSE = " << mse << std::endl;
            ofs << EbN0dB << "," << mse << std::endl;
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Mode 51 Total Execution Time: " << elapsed.count() << " seconds" << std::endl;
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
