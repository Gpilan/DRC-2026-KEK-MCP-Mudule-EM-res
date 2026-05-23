#include "TBread.h"
#include "TButility.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <chrono>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TH2.h"

#include "function.h"

namespace fs = std::filesystem;

int main(int argc, char** argv) {

    int fRunNum = std::stoi(argv[1]);
    int fMaxEvent = std::stoi(argv[2]);

    float thr = 0.3;
    float WC_calib_const = 0.05; // 0.05 mm/ns = 5 cm/microsecond

    fs::path dir("./AUX_ref");
    if (!(fs::exists(dir))) fs::create_directory(dir);

    // initialize the utility class
    TButility util = TButility();
    util.LoadMapping("../../mapping/KEK_TB2026May_PMT.root");

    // prepare CIDs that we want to use (CID = Channel ID)
    // Aux. detectors
    TBcid cid_WCX = util.GetCID("WCX"); // Wire chamber X
    TBcid cid_WCY = util.GetCID("WCY"); // Wire chamber Y
    TBcid cid_NIM = util.GetCID("NIM"); // NIM

    // prepare the histograms wa want to draw
    TH1F* hist_NIM_WCX = new TH1F("NIM_WCX", "NIM - WC X;time (ns);Evt", 2000, -800, 800);
    TH1F* hist_NIM_WCY = new TH1F("NIM_WCY", "NIM - WC Y;time (ns);Evt", 2000, -800, 800);

    TH1F* hist_NIM = new TH1F("NIM_timing", "NIM timing;time (ns);Evt", 1000, 0, 800);
    TH1F* hist_WCX = new TH1F("WCX_timing", "WC X timing;time (ns);Evt", 1000, 0, 800);
    TH1F* hist_WCY = new TH1F("WCY_timing", "WC Y timing;time (ns);Evt", 1000, 0, 800);

    TH1F* hist_pos_WCX = new TH1F("WireChamberX", "WC X position (no ref. correction);X [mm];Evt", 800, -40, 40);
    TH1F* hist_pos_WCY = new TH1F("WireChamberY", "WC Y position (no ref. correction);Y [mm];Evt", 800, -40, 40);
    TH2F* hist_WC      = new TH2F("WireChamber", "Wire Chamber (no ref. correction);X [mm];Y [mm];events", 160, -40, 40, 160, -40, 40);
    TH2F* hist_WC_fine = new TH2F("WireChamber_fine", "Wire Chamber with 0.1mm bin (no ref. correction);X [mm];Y [mm];events", 800, -40, 40, 800, -40, 40);

    // MIDs
    // 1: PMT C
    // 2: PMT S + trg
    // 8, 9, 11: MCP S
    // 11, 13, 15: MCP C
    // 16: WC, NIM
    // (Hodoscope MID 17 dropped: hodoscope not available)
    TBread<TBwaveform> readerWave = TBread<TBwaveform>(fRunNum, fMaxEvent, -1, false, "/Volumes/SSD_8TB", {16});

    // Set Maximum event
    if (fMaxEvent == -1 || fMaxEvent > readerWave.GetMaxEvent())
        fMaxEvent = readerWave.GetMaxEvent();

    for (int iEvt = 0; iEvt < fMaxEvent; iEvt++) {
        if (iEvt % 100 == 0) printProgress(iEvt, fMaxEvent);
        TBevt<TBwaveform> anEvt = readerWave.GetAnEvent();
        // Get waveform for DWCs
        std::vector<short> wave_WCX = (anEvt.GetData(cid_WCX)).waveform();
        std::vector<short> wave_WCY = (anEvt.GetData(cid_WCY)).waveform();
        std::vector<short> wave_NIM = (anEvt.GetData(cid_NIM)).waveform();

        float timing_WCX = getLeadingEdgeTime_interpolated800(wave_WCX, thr, 1, 1000);
        float timing_WCY = getLeadingEdgeTime_interpolated800(wave_WCY, thr, 1, 1000);
        float timing_NIM = getLeadingEdgeTime_interpolated800(wave_NIM, thr, 1, 1000);

        hist_WCX->Fill(timing_WCX);
        hist_WCY->Fill(timing_WCY);
        hist_NIM->Fill(timing_NIM);

        // Reference timing = NIM - WC timing
        float timeDiff_X = timing_NIM - timing_WCX;
        float timeDiff_Y = timing_NIM - timing_WCY;

        hist_NIM_WCX->Fill(timeDiff_X);
        hist_NIM_WCY->Fill(timeDiff_Y);

        // Position w/o reference correction (skip WC_X/Y_ref subtraction):
        //   pos = (NIM - WC) * calib_const,  Y is negated to match calib_DRC.cc
        float pos_X = timeDiff_X * WC_calib_const;
        float pos_Y = -1. * timeDiff_Y * WC_calib_const;

        hist_pos_WCX->Fill(pos_X);
        hist_pos_WCY->Fill(pos_Y);
        hist_WC->Fill(pos_X, pos_Y);
        hist_WC_fine->Fill(pos_X, pos_Y);
    }

    float WC_X_mean = 0.;
    float WC_Y_mean = 0.;
    float FWHM_X_ref = GetFWHM(hist_NIM_WCX, WC_X_mean); // Return FWHM of NIM - WC X, also calculates center position of X at FWHM
    float FWHM_Y_ref = GetFWHM(hist_NIM_WCY, WC_Y_mean); // Return FWHM of NIM - WC Y, also calculates center position of Y at FWHM

    std::cout << "WC X center timing at FWHM: " << WC_X_mean << std::endl;
    std::cout << "WC Y center timing at FWHM: " << WC_Y_mean << std::endl;
    std::cout << "WC X FWHM: " << FWHM_X_ref << std::endl;
    std::cout << "WC Y FWHM: " << FWHM_Y_ref << std::endl;

    std::string outName = "./AUX_ref/AUX_ref_Run_" + std::to_string(fRunNum) + ".root";
    TFile* outputRoot = new TFile(outName.c_str(), "RECREATE");
    outputRoot->cd();
    hist_NIM_WCX->Write();
    hist_NIM_WCY->Write();
    hist_NIM->Write();
    hist_WCX->Write();
    hist_WCY->Write();

    hist_pos_WCX->Write();
    hist_pos_WCY->Write();
    hist_WC->Write();
    hist_WC_fine->Write();

    outputRoot->Close();
}
