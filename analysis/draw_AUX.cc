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

namespace fs = std::filesystem;

int main(int argc, char** argv) {

    int fRunNum = std::stoi(argv[1]);
    int fMaxEvent = std::stoi(argv[2]);

    float thr = 0.3;

    int first = 135;  // Hodoscope integration range
    int last  = 270;  // Hodoscope integration range
    
    fs::path dir("./AUX_ref");   
    if (!(fs::exists(dir))) fs::create_directory(dir);

    // initialize the utility class
    TButility util = TButility();
    util.LoadMapping("../mapping/mapping_KEK.root");
 
    // prepare CIDs that we want to use (CID = Channel ID)
    // Aux. detectors
    TBcid cid_WCX = util.GetCID("WCX"); // Wire chamber X
    TBcid cid_WCY = util.GetCID("WCY"); // Wire chamber Y
    TBcid cid_NIM = util.GetCID("NIM"); // NIM

    // Hodoscope
    TBcid cid_HX1 = util.GetCID("HX1");
    TBcid cid_HX2 = util.GetCID("HX2");
    TBcid cid_HX3 = util.GetCID("HX3");
    TBcid cid_HX4 = util.GetCID("HX4");
    TBcid cid_HX5 = util.GetCID("HX5");
    TBcid cid_HX6 = util.GetCID("HX6");
    TBcid cid_HX7 = util.GetCID("HX7");
    TBcid cid_HX8 = util.GetCID("HX8");
    TBcid cid_HX9 = util.GetCID("HX9");
    TBcid cid_HX10 = util.GetCID("HX10");
    TBcid cid_HX11 = util.GetCID("HX11");
    TBcid cid_HX12 = util.GetCID("HX12");
    TBcid cid_HX13 = util.GetCID("HX13");
    TBcid cid_HX14 = util.GetCID("HX14");
    TBcid cid_HX15 = util.GetCID("HX15");
    TBcid cid_HX16 = util.GetCID("HX16");

    TBcid cid_HY1 = util.GetCID("HY1");
    TBcid cid_HY2 = util.GetCID("HY2");
    TBcid cid_HY3 = util.GetCID("HY3");
    TBcid cid_HY4 = util.GetCID("HY4");
    TBcid cid_HY5 = util.GetCID("HY5");
    TBcid cid_HY6 = util.GetCID("HY6");
    TBcid cid_HY7 = util.GetCID("HY7");
    TBcid cid_HY8 = util.GetCID("HY8");
    TBcid cid_HY9 = util.GetCID("HY9");
    TBcid cid_HY10 = util.GetCID("HY10");
    TBcid cid_HY11 = util.GetCID("HY11");
    TBcid cid_HY12 = util.GetCID("HY12");
    TBcid cid_HY13 = util.GetCID("HY13");
    TBcid cid_HY14 = util.GetCID("HY14");
    TBcid cid_HY15 = util.GetCID("HY15");
    TBcid cid_HY16 = util.GetCID("HY16");

    // prepare the histograms wa want to draw
    TH1F* hist_NIM_WCX = new TH1F("NIM_WCX", "NIM - WC X;time (ns);Evt", 2000, -800, 800);
    TH1F* hist_NIM_WCY = new TH1F("NIM_WCY", "NIM - WC Y;time (ns);Evt", 2000, -800, 800);

    TH1F* hist_NIM = new TH1F("NIM_timing", "NIM timing;time (ns);Evt", 1000, 0, 800);
    TH1F* hist_WCX = new TH1F("WCX_timing", "WC X timing;time (ns);Evt", 1000, 0, 800);
    TH1F* hist_WCY = new TH1F("WCY_timing", "WC Y timing;time (ns);Evt", 1000, 0, 800);

    TH1F* hist_Hodo_X = new TH1F("Hodo_X", "Hodo X;mm;Evt", 16, 0, 16);
    TH1F* hist_Hodo_Y = new TH1F("Hodo_Y", "Hodo Y;mm;Evt", 16, 0, 16);
    TH2F* hist_Hodoscope = new TH2F("Hodoscope", "Hodoscope;X [mm];Y [mm];events", 16, 0, 16, 16, 0, 16);

    // Change {8, 9, 13} to actual set of MIDs
    TBread<TBwaveform> readerWave = TBread<TBwaveform>(fRunNum, fMaxEvent, -1, false, "/pnfs/knu.ac.kr/data/cms/store/user/sungwon/2026_May_KEK_DRC_TB_Data", {8, 9, 13});

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

        std::vector<short> wave_HX1 = anEvt.GetData(cid_HX1).waveform();
        std::vector<short> wave_HX2 = anEvt.GetData(cid_HX2).waveform();
        std::vector<short> wave_HX3 = anEvt.GetData(cid_HX3).waveform();
        std::vector<short> wave_HX4 = anEvt.GetData(cid_HX4).waveform();
        std::vector<short> wave_HX5 = anEvt.GetData(cid_HX5).waveform();
        std::vector<short> wave_HX6 = anEvt.GetData(cid_HX6).waveform();
        std::vector<short> wave_HX7 = anEvt.GetData(cid_HX7).waveform();
        std::vector<short> wave_HX8 = anEvt.GetData(cid_HX8).waveform();
        std::vector<short> wave_HX9 = anEvt.GetData(cid_HX9).waveform();
        std::vector<short> wave_HX10 = anEvt.GetData(cid_HX10).waveform();
        std::vector<short> wave_HX11 = anEvt.GetData(cid_HX11).waveform();
        std::vector<short> wave_HX12 = anEvt.GetData(cid_HX12).waveform();
        std::vector<short> wave_HX13 = anEvt.GetData(cid_HX13).waveform();
        std::vector<short> wave_HX14 = anEvt.GetData(cid_HX14).waveform();
        std::vector<short> wave_HX15 = anEvt.GetData(cid_HX15).waveform();
        std::vector<short> wave_HX16 = anEvt.GetData(cid_HX16).waveform();

        std::vector<short> wave_HY1 = anEvt.GetData(cid_HY1).waveform();
        std::vector<short> wave_HY2 = anEvt.GetData(cid_HY2).waveform();
        std::vector<short> wave_HY3 = anEvt.GetData(cid_HY3).waveform();
        std::vector<short> wave_HY4 = anEvt.GetData(cid_HY4).waveform();
        std::vector<short> wave_HY5 = anEvt.GetData(cid_HY5).waveform();
        std::vector<short> wave_HY6 = anEvt.GetData(cid_HY6).waveform();
        std::vector<short> wave_HY7 = anEvt.GetData(cid_HY7).waveform();
        std::vector<short> wave_HY8 = anEvt.GetData(cid_HY8).waveform();
        std::vector<short> wave_HY9 = anEvt.GetData(cid_HY9).waveform();
        std::vector<short> wave_HY10 = anEvt.GetData(cid_HY10).waveform();
        std::vector<short> wave_HY11 = anEvt.GetData(cid_HY11).waveform();
        std::vector<short> wave_HY12 = anEvt.GetData(cid_HY12).waveform();
        std::vector<short> wave_HY13 = anEvt.GetData(cid_HY13).waveform();
        std::vector<short> wave_HY14 = anEvt.GetData(cid_HY14).waveform();
        std::vector<short> wave_HY15 = anEvt.GetData(cid_HY15).waveform();
        std::vector<short> wave_HY16 = anEvt.GetData(cid_HY16).waveform();

        float timing_WCX = getLeadingEdgeTime_interpolated800(wave_WCX, thr, 1, 1000);
        float timing_WCY = getLeadingEdgeTime_interpolated800(wave_WCY, thr, 1, 1000);
        float timing_NIM = getLeadingEdgeTime_interpolated800(wave_NIM, thr, 1, 1000);

        std::vector<float> intADC_HX(16);
        std::vector<float> intADC_HY(16);

        intADC_HX[0]  = GetInt(wave_HX1, first, last);
        intADC_HX[1]  = GetInt(wave_HX2, first, last);
        intADC_HX[2]  = GetInt(wave_HX3, first, last);
        intADC_HX[3]  = GetInt(wave_HX4, first, last);
        intADC_HX[4]  = GetInt(wave_HX5, first, last);
        intADC_HX[5]  = GetInt(wave_HX6, first, last);
        intADC_HX[6]  = GetInt(wave_HX7, first, last);
        intADC_HX[7]  = GetInt(wave_HX8, first, last);
        intADC_HX[8]  = GetInt(wave_HX9, first, last);
        intADC_HX[9]  = GetInt(wave_HX10, first, last);
        intADC_HX[10] = GetInt(wave_HX11, first, last);
        intADC_HX[11] = GetInt(wave_HX12, first, last);
        intADC_HX[12] = GetInt(wave_HX13, first, last);
        intADC_HX[13] = GetInt(wave_HX14, first, last);
        intADC_HX[14] = GetInt(wave_HX15, first, last);
        intADC_HX[15] = GetInt(wave_HX16, first, last);

        intADC_HY[0]  = GetInt(wave_HY1, first, last);
        intADC_HY[1]  = GetInt(wave_HY2, first, last);
        intADC_HY[2]  = GetInt(wave_HY3, first, last);
        intADC_HY[3]  = GetInt(wave_HY4, first, last);
        intADC_HY[4]  = GetInt(wave_HY5, first, last);
        intADC_HY[5]  = GetInt(wave_HY6, first, last);
        intADC_HY[6]  = GetInt(wave_HY7, first, last);
        intADC_HY[7]  = GetInt(wave_HY8, first, last);
        intADC_HY[8]  = GetInt(wave_HY9, first, last);
        intADC_HY[9]  = GetInt(wave_HY10, first, last);
        intADC_HY[10] = GetInt(wave_HY11, first, last);
        intADC_HY[11] = GetInt(wave_HY12, first, last);
        intADC_HY[12] = GetInt(wave_HY13, first, last);
        intADC_HY[13] = GetInt(wave_HY14, first, last);
        intADC_HY[14] = GetInt(wave_HY15, first, last);
        intADC_HY[15] = GetInt(wave_HY16, first, last);
        int max_X_idx_intADC = std::max_element(intADC_HX.begin(), intADC_HX.end()) - intADC_HX.begin();
        int max_Y_idx_intADC = std::max_element(intADC_HY.begin(), intADC_HY.end()) - intADC_HY.begin();

        float max_X_pos_intADC = max_X_idx_intADC + 0.5; // Position -> between 0, 16 (need to be corrected for diff calculation)
        float max_Y_pos_intADC = max_Y_idx_intADC + 0.5; // Position -> between 0, 16 (need to be corrected for diff calculation)

        hist_WCX->Fill(timing_WCX);
        hist_WCY->Fill(timing_WCY);
        hist_NIM->Fill(timing_NIM);

        // Reference timing = NIM - WC timing
        float timeDiff_X = timing_NIM - timing_WCX;
        float timeDiff_Y = timing_NIM - timing_WCY;

        hist_NIM_WCX->Fill(timeDiff_X);
        hist_NIM_WCY->Fill(timeDiff_Y);

        hist_Hodo_X->Fill(max_X_pos_intADC);
        hist_Hodo_Y->Fill(max_Y_pos_intADC);
        hist_Hodoscope->Fill(max_X_pos_intADC, max_Y_pos_intADC);
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

    hist_Hodo_X->Write();
    hist_Hodo_Y->Write();
    hist_Hodoscope->Write();

    outputRoot->Close();
}