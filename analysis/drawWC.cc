#include "TBread.h"
#include "TButility.h"

#include <filesystem>
#include <iostream>
#include <chrono>
#include <numeric>
#include <vector>
#include <fstream>
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
    float WC_calib_const = 0.05; // 0.05 mm/ns = 5 cm/microsecond

    fs::path dir("./WC");   
    if (!(fs::exists(dir))) fs::create_directory(dir);

    // initialize the utility class
    TButility util = TButility();
    util.LoadMapping("../../mapping/mapping_KEK.root");

    TFile* inputRoot = new TFile("./WC_mean/WC_mean.root");
    TH1F* hist_NIM_WCX = (TH1F*)inputRoot->Get("NIM_WCX");
    TH1F* hist_NIM_WCY = (TH1F*)inputRoot->Get("NIM_WCY");

    float WC_X_ref = 0.;
    float WC_Y_ref = 0.;
    float FWHM_X_ref = GetFWHM(hist_NIM_WCX, WC_X_ref); // Return FWHM of NIM - WC X, also calculates center position of X at FWHM
    float FWHM_Y_ref = GetFWHM(hist_NIM_WCY, WC_Y_ref); // Return FWHM of NIM - WC Y, also calculates center position of Y at FWHM

    std::cout << "WC X center timing at FWHM: " << WC_X_ref << std::endl;
    std::cout << "WC Y center timing at FWHM: " << WC_Y_ref << std::endl;
    std::cout << "WC X FWHM: " << FWHM_X_ref << std::endl;
    std::cout << "WC Y FWHM: " << FWHM_Y_ref << std::endl;
 
    // prepare CIDs that we want to use (CID = Channel ID)
    // Aux. detectors
    // Not available for now (2025-Nov-28)
    TBcid cid_WCX = util.GetCID("WCX"); // Wire chamber X
    TBcid cid_WCY = util.GetCID("WCY"); // Wire chamber Y
    TBcid cid_NIM = util.GetCID("NIM"); // NIM

    // prepare the histograms wa want to draw
    TH1F* hist_time_NIM = new TH1F("NIM_timing", "NIM timing;time (ns);Evt",  1000, 0, 800);
    TH1F* hist_time_WCX = new TH1F("WCX_timing", "WC X timing;time (ns);Evt", 1000, 0, 800);
    TH1F* hist_time_WCY = new TH1F("WCY_timing", "WC Y timing;time (ns);Evt", 1000, 0, 800);

    TH1F* hist_timeDiff_WCX = new TH1F("NIM_WCX", "NIM - WC X;time (ns);Evt", 2000, -800, 800);
    TH1F* hist_timeDiff_WCY = new TH1F("NIM_WCY", "NIM - WC Y;time (ns);Evt", 2000, -800, 800);

    TH1F* hist_timeDiff_Ref_WCX = new TH1F("Ref_NIM_WCX", "Ref X - (NIM - WC X);time (ns);Evt", 2000, -800, 800);
    TH1F* hist_timeDiff_Ref_WCY = new TH1F("Ref_NIM_WCY", "Ref Y - (NIM - WC Y);time (ns);Evt", 2000, -800, 800);

    TH1F* hist_pos_WCX = new TH1F("WireChamberX", "WC X position;X [mm];Evt", 600, -30, 30);
    TH1F* hist_pos_WCY = new TH1F("WireChamberY", "WC Y position;Y [mm];Evt", 600, -30, 30);

    TH2F* hist_WC_fine = new TH2F("WireChamber_fine", "Wire Chamber with 0.1mm bin;X [mm];Y [mm];events", 600, -30, 30, 600, -30, 30);
    TH2F* hist_WC      = new TH2F("WireChamber", "Wire Chamber;X [mm];Y [mm];events", 120, -30, 30, 120, -30, 30);
    
    // Hist after 1x1 cm^2 position cut
    TH1F* hist_time_NIM_after = new TH1F("NIM_timing_after", "NIM timing;time (ns);Evt",  1000, 0, 800);
    TH1F* hist_time_WCX_after = new TH1F("WCX_timing_after", "WC X timing;time (ns);Evt", 1000, 0, 800);
    TH1F* hist_time_WCY_after = new TH1F("WCY_timing_after", "WC Y timing;time (ns);Evt", 1000, 0, 800);

    TH1F* hist_timeDiff_WCX_after = new TH1F("NIM_WCX_after", "NIM - WC X;time (ns);Evt", 2000, -800, 800);
    TH1F* hist_timeDiff_WCY_after = new TH1F("NIM_WCY_after", "NIM - WC Y;time (ns);Evt", 2000, -800, 800);

    TH1F* hist_timeDiff_Ref_WCX_after = new TH1F("Ref_NIM_WCX_after", "Ref X - (NIM - WC X);time (ns);Evt", 2000, -800, 800);
    TH1F* hist_timeDiff_Ref_WCY_after = new TH1F("Ref_NIM_WCY_after", "Ref Y - (NIM - WC Y);time (ns);Evt", 2000, -800, 800);

    TH1F* hist_pos_WCX_after = new TH1F("WireChamberX_after", "WC X position;X [mm];Evt", 600, -30, 30);
    TH1F* hist_pos_WCY_after = new TH1F("WireChamberY_after", "WC Y position;Y [mm];Evt", 600, -30, 30);

    TH2F* hist_WC_fine_after = new TH2F("WireChamber_fine_after", "Wire Chamber with 0.1mm bin;X [mm];Y [mm];events", 600, -30, 30, 600, -30, 30);
    TH2F* hist_WC_after      = new TH2F("WireChamber_after", "Wire Chamber;X [mm];Y [mm];events", 120, -30, 30, 120, -30, 30);

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

        // Reference timing = NIM - WC timing
        float timeDiff_NIM_X = timing_NIM - timing_WCX;
        float timeDiff_NIM_Y = timing_NIM - timing_WCY;

        // We do not know which timing corresponds to which side (L or R, U or D)
        // By applying X cut [0, 5] mm -> found that Diff > 0 = Right
        // By applying Y cut [0, 5] mm -> found that Diff > 0 = Down
        float timeDiff_Ref_X = WC_X_ref - timeDiff_NIM_X; // diff > 0 -> right side, diff < 0 -> left side
        float timeDiff_Ref_Y = WC_Y_ref - timeDiff_NIM_Y; // diff > 0 -> down side, diff < 0 -> up side

        float pos_X = timeDiff_Ref_X * WC_calib_const; // diff > 0 -> right side, diff < 0 -> left side
        float pos_Y = -1. * timeDiff_Ref_Y * WC_calib_const; // diff > 0 -> down side, diff < 0 -> up side

        hist_time_WCX->Fill(timing_WCX);
        hist_time_WCY->Fill(timing_WCY);
        hist_time_NIM->Fill(timing_NIM);

        hist_timeDiff_WCX->Fill(timeDiff_NIM_X);
        hist_timeDiff_WCY->Fill(timeDiff_NIM_Y);

        hist_timeDiff_Ref_WCX->Fill(timeDiff_Ref_X);
        hist_timeDiff_Ref_WCY->Fill(timeDiff_Ref_Y);

        hist_pos_WCX->Fill(pos_X);
        hist_pos_WCY->Fill(pos_Y);
        hist_WC->Fill(pos_X, pos_Y);
        hist_WC_fine->Fill(pos_X, pos_Y);
    }


    std::string outFile = "./WC/WC_Run_" + std::to_string(fRunNum) + ".root";
    TFile* outputRoot = new TFile(outFile.c_str(), "RECREATE");
    outputRoot->cd();

    hist_time_WCX->Write();
    hist_time_WCY->Write();
    hist_time_NIM->Write();
    hist_timeDiff_WCX->Write();
    hist_timeDiff_WCY->Write();
    hist_timeDiff_Ref_WCX->Write();
    hist_timeDiff_Ref_WCY->Write();
    hist_pos_WCX->Write();
    hist_pos_WCY->Write();
    hist_WC->Write();
    hist_WC_fine->Write();

    outputRoot->Close();
}