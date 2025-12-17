#include "TBread.h"
#include "TButility.h"

#include <filesystem>
#include <iostream>
#include <chrono>
#include <numeric>
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

namespace fs = std::__fs::filesystem;

int main(int argc, char** argv) {

    int fRunNum = std::stoi(argv[1]);
    int fMaxEvent = std::stoi(argv[2]);
    int fMaxFile = -1;

    float thr = 0.3;

    fs::path dir("./WC_mean");   
    if (!(fs::exists(dir))) fs::create_directory(dir);

    // initialize the utility class
    TButility util = TButility();
    util.LoadMapping("../../mapping/mapping_KEK.root");
 
    // prepare CIDs that we want to use (CID = Channel ID)
    // Aux. detectors
    // Not available for now (2025-Nov-28)
    TBcid cid_WCX = util.GetCID("WCX"); // Wire chamber X
    TBcid cid_WCY = util.GetCID("WCY"); // Wire chamber Y
    TBcid cid_NIM = util.GetCID("NIM"); // NIM

    // prepare the histograms wa want to draw
    TH1F* hist_NIM_WCX = new TH1F("NIM_WCX", "NIM - WC X;time (ns);Evt", 2000, -800, 800);
    TH1F* hist_NIM_WCY = new TH1F("NIM_WCY", "NIM - WC Y;time (ns);Evt", 2000, -800, 800);

    TH1F* hist_NIM = new TH1F("NIM_timing", "NIM timing;time (ns);Evt", 1000, 0, 800);
    TH1F* hist_WCX = new TH1F("WCX_timing", "WC X timing;time (ns);Evt", 1000, 0, 800);
    TH1F* hist_WCY = new TH1F("WCY_timing", "WC Y timing;time (ns);Evt", 1000, 0, 800);

    TBread<TBwaveform> readerWave = TBread<TBwaveform>(fRunNum, fMaxEvent, fMaxFile, false, "/Volumes/SSD_8TB", {13});

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
    }

    float WC_X_mean = 0.;
    float WC_Y_mean = 0.;
    float FWHM_X_ref = GetFWHM(hist_NIM_WCX, WC_X_mean); // Return FWHM of NIM - WC X, also calculates center position of X at FWHM
    float FWHM_Y_ref = GetFWHM(hist_NIM_WCY, WC_Y_mean); // Return FWHM of NIM - WC Y, also calculates center position of Y at FWHM

    std::cout << "WC X center timing at FWHM: " << WC_X_mean << std::endl;
    std::cout << "WC Y center timing at FWHM: " << WC_Y_mean << std::endl;
    std::cout << "WC X FWHM: " << FWHM_X_ref << std::endl;
    std::cout << "WC Y FWHM: " << FWHM_Y_ref << std::endl;

    std::string outName = "./WC_mean/WC_mean.root";
    TFile* outputRoot = new TFile(outName.c_str(), "RECREATE");
    outputRoot->cd();
    hist_NIM_WCX->Write();
    hist_NIM_WCY->Write();
    hist_NIM->Write();
    hist_WCX->Write();
    hist_WCY->Write();

    outputRoot->Close();
}