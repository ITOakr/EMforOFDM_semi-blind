#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include "simulator_misc.h"
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
double avgPower;
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

    SimulatorMisc sim(params);

    int mode_select = 0;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Select Miscellaneous Simulation / Utility Mode" << std::endl;
    std::cout << "4: Average power of the true channel response H" << std::endl;
    std::cout << "10: CSV output of channel magnitude response |H(k, l)|" << std::endl;
    std::cout << "23: Export transmit waveform X over time (k=10)" << std::endl;
    std::cout << "24: Export faded waveform HX over time (k=10)" << std::endl;
    std::cout << "25: Export channel magnitude |H| over time (k=10)" << std::endl;
    std::cout << "26: Export frequency response |H(k, l)| along k for fixed l" << std::endl;
    std::cout << "27: Export average impulse response along q for fixed l" << std::endl;
    std::cout << "28: Export estimated impulse response (l=0, Q=16) to CSV" << std::endl;
    std::cout << "30: Sweep CRLB-based MSE vs Eb/N0 (fixed Doppler)" << std::endl;
    std::cout << "31: Simulate impulse-response h MSE vs Eb/N0" << std::endl;
    std::cout << "32: Known-model, known-noise ML validation (H MSE vs Eb/N0)" << std::endl;
    std::cout << "33: Known-model, known-noise ML validation (H MSE at l=0 vs Eb/N0)" << std::endl;
    std::cout << "34: Instantaneous SNR (gamma) mean vs Eb/N0 at (l,k)=(0,0)" << std::endl;
    std::cout << "35: Instantaneous SNR (gamma) vs Eb/N0 at (l,k)=(0,0)" << std::endl;
    std::cout << "36: Instantaneous SNR (gamma) vs Doppler at (l,k)=(0,0)" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> mode_select;

    std::cout << "Enter number of trials:";
    std::cin >> numberOfTrials;
    sim.setTrialNum(numberOfTrials);

    std::string modulationName = getModulationSchemeName(params.NUMBER_OF_BIT);
    std::string timeStr = getCurrentTimeString();

    if (mode_select == 4)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        sim.setDopplerFrequency(dopplerFrequency);
        
        avgPower = sim.getAveragePower_simulation();
        
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Average power of the true channel response H over " << numberOfTrials << " trials." << std::endl;
        std::cout << "f_d*T_s = " << dopplerFrequency << ", Average Power = " << avgPower << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        return 0; 
    }
    else if (mode_select == 10)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s for single trial:" ;
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_Channel_Magnitude_Response.csv";
        ofs.open(fileName);

        sim.setDopplerFrequency(dopplerFrequency);
        sim.saveChannelMagnitudeResponseToCSV(ofs, dopplerFrequency);
        
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Channel magnitude response saved to: " << fileName << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
    }
    else if (mode_select == 23)
    {
        fileName = outputDir + timeStr + "_" + modulationName + "_TxWaveform_k10.csv";
        int target_k = 10;
        
        std::cout << "Exporting Tx Signal for k=" << target_k << " ..." << std::endl;
        sim.runExportTxWaveform(target_k, fileName);
    }
    else if (mode_select == 24)
    {
        double inputDoppler;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> inputDoppler;
        sim.setDopplerFrequency(inputDoppler);

        fileName = outputDir + timeStr + "_" + modulationName + "_FadedWaveform_k10.csv";
        int target_k = 10;

        std::cout << "Exporting Faded Signal for k=" << target_k << " ..." << std::endl;
        sim.runExportFadedWaveform(target_k, fileName);
    }
    else if (mode_select == 25)
    {
        double inputDoppler;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> inputDoppler;
        sim.setDopplerFrequency(inputDoppler);

        int target_k;
        std::cout << "Enter target subcarrier index (k): ";
        std::cin >> target_k;

        fileName = outputDir + timeStr + "_" + modulationName + "_ChannelMagnitude_k" + std::to_string(target_k) + ".csv";
        
        std::cout << "Exporting Channel Magnitude for k=" << target_k << " ..." << std::endl;
        sim.runExportChannelMagnitude(target_k, fileName);
    }
    else if (mode_select == 26)
    {
        double inputDoppler;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> inputDoppler;
        sim.setDopplerFrequency(inputDoppler);

        int target_l;
        std::cout << "Enter target symbol index (l): ";
        std::cin >> target_l;

        fileName = outputDir + timeStr + "_" + modulationName + "_FreqResp_k_at_l" + std::to_string(target_l) + ".csv";
        
        std::cout << "Exporting Frequency Response vs k..." << std::endl;
        sim.saveFrequencyResponseByK(target_l, fileName);
        std::cout << "Saved to: " << fileName << std::endl;
    }
    else if (mode_select == 27)
    {
        double inputDoppler;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> inputDoppler;
        sim.setDopplerFrequency(inputDoppler);

        int target_l = 0;

        fileName = outputDir + timeStr + "_" + modulationName + "_AverageImpulseResp_q.csv";
        sim.saveAverageImpulseResponseByQ(target_l, fileName);
        
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Average impulse response saved to: " << fileName << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
    }
    else if (mode_select == 28)
    {
        double inputDoppler;
        std::cout << "Enter normalized Doppler frequency (f_d T_s):";
        std::cin >> inputDoppler;
        sim.setDopplerFrequency(inputDoppler);

        int fixedEbN0dB;
        std::cout << "Enter fixed Eb/N0 [dB]:";
        std::cin >> fixedEbN0dB;
        sim.setNoiseSD(fixedEbN0dB);

        fileName = outputDir + timeStr + "_" + modulationName + "_EbN0_" + std::to_string(fixedEbN0dB) + "_EstImpulseResponse_l0.csv";
        ofs.open(fileName);

        std::cout << "Exporting Estimated Impulse Response (l=0, Q=16)..." << std::endl;
        sim.saveEstimatedImpulseResponseToCSV(ofs, inputDoppler);
        std::cout << "Saved to: " << fileName << std::endl;
    }
    else if (mode_select == 30)
    {
        fileName = outputDir + timeStr + "_" + modulationName + "_CRLB_vs_EbN0.csv";
        ofs.open(fileName);
        ofs << "EbN0dB,CRLB_MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setNoiseSD(EbN0dB);
            double crlb_mse = sim.getTheoreticalCRLB_H_MSE_FinalForm();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", CRLB MSE = " << crlb_mse << std::endl;
            ofs << EbN0dB << "," << crlb_mse << std::endl;
        }
    }
    else if (mode_select == 31)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "_fdTs_" + std::to_string(dopplerFrequency) + "_ImpulseResponse_MSE_vs_EbN0.csv";
        ofs.open(fileName);
        ofs << "EbN0dB,Impulse_MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            mse = sim.getImpulseResponseMSE_simulation();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", Impulse MSE = " << mse << std::endl;
            ofs << EbN0dB << "," << mse << std::endl;
        }
    }
    else if (mode_select == 32)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_KnownModelKnownNoise_HMSE_vs_EbN0.csv";
        ofs.open(fileName);
        ofs << "EbN0dB,H_MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double H_mse = sim.get_H_MSE_with_known_model_and_noise_during_pilot();

            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", H_MSE = " << H_mse << std::endl;
            ofs << EbN0dB << "," << H_mse << std::endl;
        }
    }
    else if (mode_select == 33)
    {
        double dopplerFrequency;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_KnownModelKnownNoise_HMSE_vs_EbN0.csv";
        ofs.open(fileName);
        ofs << "EbN0dB,H_MSE" << std::endl;

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);

            double H_mse = sim.get_H_MSE_with_known_model_and_noise_at_l0();

            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", H_MSE = " << H_mse << std::endl;
            ofs << EbN0dB << "," << H_mse << std::endl;
        }
    }
    else if (mode_select == 34)
    {
        double dopplerFrequency;
        double gamma_mean;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_Gamma_vs_EbN0_MODE34.csv";
        ofs.open(fileName);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            gamma_mean = sim.getInstantaneousSNR_Mode34_simulation();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", Mean Gamma(l=0,k=0) = " << gamma_mean << std::endl;
            ofs << EbN0dB << "," << gamma_mean << std::endl;
        }
    }
    else if (mode_select == 35)
    {
        double dopplerFrequency;
        double gamma_mean;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_Gamma_vs_EbN0_MODE35.csv";
        ofs.open(fileName);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            gamma_mean = sim.getInstantaneousSNR_Mode35_simulation();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", Mean Gamma(l=0,k=0) = " << gamma_mean << std::endl;
            ofs << EbN0dB << "," << gamma_mean << std::endl;
        }
    }
    else if (mode_select == 36)
    {
        double dopplerFrequency;
        double gamma_mean;
        std::cout << "Enter normalized Doppler f_d*T_s:";
        std::cin >> dopplerFrequency;

        fileName = outputDir + timeStr + "_" + modulationName + "f_dT_s =" + std::to_string(dopplerFrequency) + "_Gamma_vs_EbN0_MODE36.csv";
        ofs.open(fileName);

        for (int EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
            sim.setDopplerFrequency(dopplerFrequency);
            sim.setNoiseSD(EbN0dB);
            
            gamma_mean = sim.getInstantaneousSNR_Mode36_simulation();
            
            std::cout << "-----------" << std::endl;
            std::cout << "EbN0dB = " << EbN0dB << ", Mean Gamma(l=0,k=0) = " << gamma_mean << std::endl;
            ofs << EbN0dB << "," << gamma_mean << std::endl;
        }
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
