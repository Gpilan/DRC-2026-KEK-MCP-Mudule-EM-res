#include "TBread.h"
#include "TButility.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TF1.h"

#include "function.h"

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {

    // setup for prompt analysis
    // it could be like
    int fRunNum = std::stoi(argv[1]);
    int fMaxEvent = std::stoi(argv[2]);

    fs::path dir("./Calib");
    if (!(fs::exists(dir))) fs::create_directory(dir);

    // initialize the utility class
    TButility util = TButility();
    util.LoadMapping("../../mapping/KEK_TB2026May_PMT.root");

    int C_first = 140; // Module integration range
    int C_last = 280;  // Module integration range
    int S_first = 145; // Module integration range
    int S_last = 340;  // Module integration range

    int C_first_T5 = 150.;   // integration range for T5
    int C_last_T5 = 290.; // integration range for T5
    int S_first_T5 = 160.;   // integration range for T5
    int S_last_T5 = 360.; // integration range for T5

    float cut_WC_X_min = -5; // mm, trigger has size of 1x1 cm^2
    float cut_WC_X_max = 5;  // mm, trigger has size of 1x1 cm^2
    float cut_WC_Y_min = -5; // mm, trigger has size of 1x1 cm^2
    float cut_WC_Y_max = 5;  // mm, trigger has size of 1x1 cm^2

    float leadingEdge_thr = 0.3;
    float WC_calib_const = 0.05; // 0.05 mm/ns = 5 cm/microsecond
    TFile* AUX_rootfile = new TFile(("./AUX_ref/AUX_ref_Run_" + std::to_string(fRunNum) + ".root").c_str());
    TH1F* hist_NIM_WCX = (TH1F*)AUX_rootfile->Get("NIM_WCX");
    TH1F* hist_NIM_WCY = (TH1F*)AUX_rootfile->Get("NIM_WCY");
    float WC_X_ref = 0.;
    float WC_Y_ref = 0.;
    GetFWHM(hist_NIM_WCX, WC_X_ref); // Return FWHM of NIM - WC X, also calculates center position of X at FWHM
    GetFWHM(hist_NIM_WCY, WC_Y_ref); // Return FWHM of NIM - WC Y, also calculates center position of Y at FWHM
    std::cout << "WC X center timing at FWHM: " << WC_X_ref << std::endl;
    std::cout << "WC Y center timing at FWHM: " << WC_Y_ref << std::endl;

    // Scaled factor for each channel
    float scaleFactor_C = 0.9846;
    float scaleFactor_S = 1.0025;

    float calib_T1C = 0.000063674;
    float calib_T2C = 0.000064320;
    float calib_T3C = 0.000064528;
    float calib_T4C = 0.000066421;
    float calib_T5C = 0.000060354;
    float calib_T6C = 0.000065821;
    float calib_T7C = 0.000064138;
    float calib_T8C = 0.000060189;
    float calib_T9C = 0.000061569;

    float calib_T1S = 0.000028695;
    float calib_T2S = 0.000030289;
    float calib_T3S = 0.000028662;
    float calib_T4S = 0.000028088;
    float calib_T5S = 0.000028435;
    float calib_T6S = 0.000028842;
    float calib_T7S = 0.000027692;
    float calib_T8S = 0.000028358;
    float calib_T9S = 0.000029474;

    // prepare C    // prepare CIDs that we want to use (CID = Channel ID)
    // 3x3 tower CIDs
    TBcid cid_T1C = util.GetCID("T1-C");
    TBcid cid_T2C = util.GetCID("T2-C");
    TBcid cid_T3C = util.GetCID("T3-C");
    TBcid cid_T4C = util.GetCID("T4-C");
    TBcid cid_T5C = util.GetCID("T5-C");
    TBcid cid_T6C = util.GetCID("T6-C");
    TBcid cid_T7C = util.GetCID("T7-C");
    TBcid cid_T8C = util.GetCID("T8-C");
    TBcid cid_T9C = util.GetCID("T9-C");

    TBcid cid_T1S = util.GetCID("T1-S");
    TBcid cid_T2S = util.GetCID("T2-S");
    TBcid cid_T3S = util.GetCID("T3-S");
    TBcid cid_T4S = util.GetCID("T4-S");
    TBcid cid_T5S = util.GetCID("T5-S");
    TBcid cid_T6S = util.GetCID("T6-S");
    TBcid cid_T7S = util.GetCID("T7-S");
    TBcid cid_T8S = util.GetCID("T8-S");
    TBcid cid_T9S = util.GetCID("T9-S");
    // Aux. detectors
    TBcid cid_WCX = util.GetCID("WCX"); // Wire chamber X
    TBcid cid_WCY = util.GetCID("WCY"); // Wire chamber Y
    TBcid cid_NIM = util.GetCID("NIM"); // NIM

    // prepare the histograms wa want to draw
    TH1F* hist_T1C_intADC = new TH1F("T1C_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T2C_intADC = new TH1F("T2C_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T3C_intADC = new TH1F("T3C_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T4C_intADC = new TH1F("T4C_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T5C_intADC = new TH1F("T5C_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T6C_intADC = new TH1F("T6C_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T7C_intADC = new TH1F("T7C_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T8C_intADC = new TH1F("T8C_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T9C_intADC = new TH1F("T9C_intADC" , ";intADC;nEvents", 320, -20000, 300000);

    TH1F* hist_T1S_intADC = new TH1F("T1S_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T2S_intADC = new TH1F("T2S_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T3S_intADC = new TH1F("T3S_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T4S_intADC = new TH1F("T4S_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T5S_intADC = new TH1F("T5S_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T6S_intADC = new TH1F("T6S_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T7S_intADC = new TH1F("T7S_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T8S_intADC = new TH1F("T8S_intADC" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T9S_intADC = new TH1F("T9S_intADC" , ";intADC;nEvents", 320, -20000, 300000);

    TH1F* hist_T1C_peakADC = new TH1F("T1C_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T2C_peakADC = new TH1F("T2C_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T3C_peakADC = new TH1F("T3C_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T4C_peakADC = new TH1F("T4C_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T5C_peakADC = new TH1F("T5C_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T6C_peakADC = new TH1F("T6C_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T7C_peakADC = new TH1F("T7C_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T8C_peakADC = new TH1F("T8C_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T9C_peakADC = new TH1F("T9C_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);

    TH1F* hist_T1S_peakADC = new TH1F("T1S_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T2S_peakADC = new TH1F("T2S_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T3S_peakADC = new TH1F("T3S_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T4S_peakADC = new TH1F("T4S_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T5S_peakADC = new TH1F("T5S_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T6S_peakADC = new TH1F("T6S_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T7S_peakADC = new TH1F("T7S_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T8S_peakADC = new TH1F("T8S_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T9S_peakADC = new TH1F("T9S_peakADC" , ";peakADC;nEvents", 4096, 0, 4096);

    TH1F* hist_Edep_T1C = new TH1F("Edep_T1C" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T2C = new TH1F("Edep_T2C" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T3C = new TH1F("Edep_T3C" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T4C = new TH1F("Edep_T4C" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T5C = new TH1F("Edep_T5C" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T6C = new TH1F("Edep_T6C" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T7C = new TH1F("Edep_T7C" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T8C = new TH1F("Edep_T8C" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T9C = new TH1F("Edep_T9C" , ";GeV;Event", 110, -1, 10);

    TH1F* hist_Edep_T1S = new TH1F("Edep_T1S" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T2S = new TH1F("Edep_T2S" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T3S = new TH1F("Edep_T3S" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T4S = new TH1F("Edep_T4S" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T5S = new TH1F("Edep_T5S" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T6S = new TH1F("Edep_T6S" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T7S = new TH1F("Edep_T7S" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T8S = new TH1F("Edep_T8S" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T9S = new TH1F("Edep_T9S" , ";GeV;Event", 110, -1, 10);

    TH1F* hist_totalEdep_C          = new TH1F("totalEdep_C",    ";GeV;Event", 110, -1, 10);
    TH1F* hist_totalEdep_S          = new TH1F("totalEdep_S",    ";GeV;Event", 110, -1, 10);
    TH1F* hist_totalEdep_Comb       = new TH1F("totalEdep_Comb", ";GeV;Event", 110, -2, 20);

    TH1F* hist_totalEdepScaled_C    = new TH1F("totalEdepScaled_C",    ";GeV;Event", 110, -1, 10);
    TH1F* hist_totalEdepScaled_S    = new TH1F("totalEdepScaled_S",    ";GeV;Event", 110, -1, 10);
    TH1F* hist_totalEdepScaled_Comb = new TH1F("totalEdepScaled_Comb", ";GeV;Event", 110, -2, 20);

    // Wire chamber 1D position plot
    TH1F* hist_pos_WCX = new TH1F("WireChamberX", "WC X position;X [mm];Evt", 600, -30, 30);
    TH1F* hist_pos_WCY = new TH1F("WireChamberY", "WC Y position;Y [mm];Evt", 600, -30, 30);

    // Wire chamber 2D position plot
    TH2F* hist_WC_fine = new TH2F("WireChamber_fine", "Wire Chamber with 0.1mm bin;X [mm];Y [mm];events", 600, -30, 30, 600, -30, 30);
    TH2F* hist_WC      = new TH2F("WireChamber", "Wire Chamber;X [mm];Y [mm];events", 120, -30, 30, 120, -30, 30);

    // Histogram after event selection
    TH1F* hist_T1C_intADC_after = new TH1F("T1C_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T2C_intADC_after = new TH1F("T2C_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T3C_intADC_after = new TH1F("T3C_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T4C_intADC_after = new TH1F("T4C_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T5C_intADC_after = new TH1F("T5C_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T6C_intADC_after = new TH1F("T6C_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T7C_intADC_after = new TH1F("T7C_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T8C_intADC_after = new TH1F("T8C_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T9C_intADC_after = new TH1F("T9C_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);

    TH1F* hist_T1S_intADC_after = new TH1F("T1S_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T2S_intADC_after = new TH1F("T2S_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T3S_intADC_after = new TH1F("T3S_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T4S_intADC_after = new TH1F("T4S_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T5S_intADC_after = new TH1F("T5S_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T6S_intADC_after = new TH1F("T6S_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T7S_intADC_after = new TH1F("T7S_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T8S_intADC_after = new TH1F("T8S_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);
    TH1F* hist_T9S_intADC_after = new TH1F("T9S_intADC_after" , ";intADC;nEvents", 320, -20000, 300000);

    TH1F* hist_T1C_peakADC_after = new TH1F("T1C_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T2C_peakADC_after = new TH1F("T2C_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T3C_peakADC_after = new TH1F("T3C_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T4C_peakADC_after = new TH1F("T4C_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T5C_peakADC_after = new TH1F("T5C_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T6C_peakADC_after = new TH1F("T6C_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T7C_peakADC_after = new TH1F("T7C_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T8C_peakADC_after = new TH1F("T8C_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T9C_peakADC_after = new TH1F("T9C_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);

    TH1F* hist_T1S_peakADC_after = new TH1F("T1S_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T2S_peakADC_after = new TH1F("T2S_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T3S_peakADC_after = new TH1F("T3S_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T4S_peakADC_after = new TH1F("T4S_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T5S_peakADC_after = new TH1F("T5S_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T6S_peakADC_after = new TH1F("T6S_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T7S_peakADC_after = new TH1F("T7S_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T8S_peakADC_after = new TH1F("T8S_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);
    TH1F* hist_T9S_peakADC_after = new TH1F("T9S_peakADC_after" , ";peakADC;nEvents", 4096, 0, 4096);

    TH1F* hist_Edep_T1C_after = new TH1F("Edep_T1C_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T2C_after = new TH1F("Edep_T2C_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T3C_after = new TH1F("Edep_T3C_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T4C_after = new TH1F("Edep_T4C_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T5C_after = new TH1F("Edep_T5C_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T6C_after = new TH1F("Edep_T6C_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T7C_after = new TH1F("Edep_T7C_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T8C_after = new TH1F("Edep_T8C_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T9C_after = new TH1F("Edep_T9C_after" , ";GeV;Event", 110, -1, 10);

    TH1F* hist_Edep_T1S_after = new TH1F("Edep_T1S_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T2S_after = new TH1F("Edep_T2S_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T3S_after = new TH1F("Edep_T3S_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T4S_after = new TH1F("Edep_T4S_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T5S_after = new TH1F("Edep_T5S_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T6S_after = new TH1F("Edep_T6S_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T7S_after = new TH1F("Edep_T7S_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T8S_after = new TH1F("Edep_T8S_after" , ";GeV;Event", 110, -1, 10);
    TH1F* hist_Edep_T9S_after = new TH1F("Edep_T9S_after" , ";GeV;Event", 110, -1, 10);

    TH1F* hist_totalEdep_C_after          = new TH1F("totalEdep_C_after",    ";GeV;Event", 110, -1, 10);
    TH1F* hist_totalEdep_S_after          = new TH1F("totalEdep_S_after",    ";GeV;Event", 110, -1, 10);
    TH1F* hist_totalEdep_Comb_after       = new TH1F("totalEdep_Comb_after", ";GeV;Event", 110, -2, 20);

    TH1F* hist_totalEdepScaled_C_after    = new TH1F("totalEdepScaled_C_after",    ";GeV;Event", 110, -1, 10);
    TH1F* hist_totalEdepScaled_S_after    = new TH1F("totalEdepScaled_S_after",    ";GeV;Event", 110, -1, 10);
    TH1F* hist_totalEdepScaled_Comb_after = new TH1F("totalEdepScaled_Comb_after", ";GeV;Event", 110, -2, 20);

    TH1F* hist_pos_WCX_after = new TH1F("WireChamberX_after", "WC X position;X [mm];Evt", 600, -30, 30);
    TH1F* hist_pos_WCY_after = new TH1F("WireChamberY_after", "WC Y position;Y [mm];Evt", 600, -30, 30);

    TH2F* hist_WC_fine_after = new TH2F("WireChamber_fine_after", "Wire Chamber with 0.1mm bin;X [mm];Y [mm];events", 600, -30, 30, 600, -30, 30);
    TH2F* hist_WC_after      = new TH2F("WireChamber_after", "Wire Chamber;X [mm];Y [mm];events", 120, -30, 30, 120, -30, 30);

    // MIDs
    // 1: PMT C
    // 2: PMT S + trg
    // 8, 9, 11: MCP S
    // 11, 13, 15: MCP C
    // 16: WC, NIM
    // (Hodoscope MID 17 dropped: hodoscope not available)
    TBread<TBwaveform> readerWave = TBread<TBwaveform>(fRunNum, fMaxEvent, -1, false, "/Volumes/SSD_8TB", {1, 2, 16});

    // Set Maximum event
    if (fMaxEvent == -1 || fMaxEvent > readerWave.GetMaxEvent())
        fMaxEvent = readerWave.GetMaxEvent();

    for (int iEvt = 0; iEvt < fMaxEvent; iEvt++) {
        if (iEvt % 100 == 0) printProgress(iEvt, fMaxEvent);

        // Load event
        TBevt<TBwaveform> anEvt = readerWave.GetAnEvent();

        //////////////////////////////////////////////////////////////////////
        // Preparing the waveform for each channel
        //////////////////////////////////////////////////////////////////////
        std::vector<short> wave_T1C = anEvt.GetData(cid_T1C).waveform();
        std::vector<short> wave_T2C = anEvt.GetData(cid_T2C).waveform();
        std::vector<short> wave_T3C = anEvt.GetData(cid_T3C).waveform();
        std::vector<short> wave_T4C = anEvt.GetData(cid_T4C).waveform();
        std::vector<short> wave_T5C = anEvt.GetData(cid_T5C).waveform();
        std::vector<short> wave_T6C = anEvt.GetData(cid_T6C).waveform();
        std::vector<short> wave_T7C = anEvt.GetData(cid_T7C).waveform();
        std::vector<short> wave_T8C = anEvt.GetData(cid_T8C).waveform();
        std::vector<short> wave_T9C = anEvt.GetData(cid_T9C).waveform();

        std::vector<short> wave_T1S = anEvt.GetData(cid_T1S).waveform();
        std::vector<short> wave_T2S = anEvt.GetData(cid_T2S).waveform();
        std::vector<short> wave_T3S = anEvt.GetData(cid_T3S).waveform();
        std::vector<short> wave_T4S = anEvt.GetData(cid_T4S).waveform();
        std::vector<short> wave_T5S = anEvt.GetData(cid_T5S).waveform();
        std::vector<short> wave_T6S = anEvt.GetData(cid_T6S).waveform();
        std::vector<short> wave_T7S = anEvt.GetData(cid_T7S).waveform();
        std::vector<short> wave_T8S = anEvt.GetData(cid_T8S).waveform();
        std::vector<short> wave_T9S = anEvt.GetData(cid_T9S).waveform();

        std::vector<short> wave_WCX = (anEvt.GetData(cid_WCX)).waveform();
        std::vector<short> wave_WCY = (anEvt.GetData(cid_WCY)).waveform();

        std::vector<short> wave_NIM = (anEvt.GetData(cid_NIM)).waveform();

        //////////////////////////////////////////////////////////////////////
        // Position measurement - Wirechamber
        //////////////////////////////////////////////////////////////////////
        // Get waveform for Wire chambers and NIM
        float timing_WCX = getLeadingEdgeTime_interpolated800(wave_WCX, leadingEdge_thr, 1, 1000);
        float timing_WCY = getLeadingEdgeTime_interpolated800(wave_WCY, leadingEdge_thr, 1, 1000);
        float timing_NIM = getLeadingEdgeTime_interpolated800(wave_NIM, leadingEdge_thr, 1, 1000);
        // Reference timing = NIM - WC timing
        float timeDiff_NIM_X = timing_NIM - timing_WCX;
        float timeDiff_NIM_Y = timing_NIM - timing_WCY;
        // We do not know which timing corresponds to which side (L or R, U or D)
        float timeDiff_Ref_X = WC_X_ref - timeDiff_NIM_X; // diff > 0 -> right side, diff < 0 -> left side
        float timeDiff_Ref_Y = WC_Y_ref - timeDiff_NIM_Y; // diff > 0 -> down side, diff < 0 -> up side
        // Position = (Time difference btw reference and (NIM-MC)) * WC calib. const.
        float pos_X = timeDiff_Ref_X * WC_calib_const;       // diff > 0 -> right side, diff < 0 -> left side
        float pos_Y = -1. * timeDiff_Ref_Y * WC_calib_const; // diff > 0 -> down side, diff < 0 -> up side

        //////////////////////////////////////////////////////////////////////
        // IntADC and Energy measurement
        //////////////////////////////////////////////////////////////////////
        float intADC_T1C = GetInt(wave_T1C, C_first, C_last);
        float intADC_T2C = GetInt(wave_T2C, C_first, C_last);
        float intADC_T3C = GetInt(wave_T3C, C_first, C_last);
        float intADC_T4C = GetInt(wave_T4C, C_first, C_last);
        float intADC_T5C = GetInt(wave_T5C, C_first_T5, C_last_T5);
        float intADC_T6C = GetInt(wave_T6C, C_first, C_last);
        float intADC_T7C = GetInt(wave_T7C, C_first, C_last);
        float intADC_T8C = GetInt(wave_T8C, C_first, C_last);
        float intADC_T9C = GetInt(wave_T9C, C_first, C_last);

        float intADC_T1S = GetInt(wave_T1S, S_first, S_last);
        float intADC_T2S = GetInt(wave_T2S, S_first, S_last);
        float intADC_T3S = GetInt(wave_T3S, S_first, S_last);
        float intADC_T4S = GetInt(wave_T4S, S_first, S_last);
        float intADC_T5S = GetInt(wave_T5S, S_first_T5, S_last_T5);
        float intADC_T6S = GetInt(wave_T6S, S_first, S_last);
        float intADC_T7S = GetInt(wave_T7S, S_first, S_last);
        float intADC_T8S = GetInt(wave_T8S, S_first, S_last);
        float intADC_T9S = GetInt(wave_T9S, S_first, S_last);

        float peakADC_T1C = GetPeak(wave_T1C, 100, 700);
        float peakADC_T2C = GetPeak(wave_T2C, 100, 700);
        float peakADC_T3C = GetPeak(wave_T3C, 100, 700);
        float peakADC_T4C = GetPeak(wave_T4C, 100, 700);
        float peakADC_T5C = GetPeak(wave_T5C, 100, 700);
        float peakADC_T6C = GetPeak(wave_T6C, 100, 700);
        float peakADC_T7C = GetPeak(wave_T7C, 100, 700);
        float peakADC_T8C = GetPeak(wave_T8C, 100, 700);
        float peakADC_T9C = GetPeak(wave_T9C, 100, 700);

        float peakADC_T1S = GetPeak(wave_T1S, 100, 700);
        float peakADC_T2S = GetPeak(wave_T2S, 100, 700);
        float peakADC_T3S = GetPeak(wave_T3S, 100, 700);
        float peakADC_T4S = GetPeak(wave_T4S, 100, 700);
        float peakADC_T5S = GetPeak(wave_T5S, 100, 700);
        float peakADC_T6S = GetPeak(wave_T6S, 100, 700);
        float peakADC_T7S = GetPeak(wave_T7S, 100, 700);
        float peakADC_T8S = GetPeak(wave_T8S, 100, 700);
        float peakADC_T9S = GetPeak(wave_T9S, 100, 700);

        float Edep_T1C = (intADC_T1C * calib_T1C);
        float Edep_T2C = (intADC_T2C * calib_T2C);
        float Edep_T3C = (intADC_T3C * calib_T3C);
        float Edep_T4C = (intADC_T4C * calib_T4C);
        float Edep_T5C = (intADC_T5C * calib_T5C);
        float Edep_T6C = (intADC_T6C * calib_T6C);
        float Edep_T7C = (intADC_T7C * calib_T7C);
        float Edep_T8C = (intADC_T8C * calib_T8C);
        float Edep_T9C = (intADC_T9C * calib_T9C);

        float Edep_T1S = (intADC_T1S * calib_T1S);
        float Edep_T2S = (intADC_T2S * calib_T2S);
        float Edep_T3S = (intADC_T3S * calib_T3S);
        float Edep_T4S = (intADC_T4S * calib_T4S);
        float Edep_T5S = (intADC_T5S * calib_T5S);
        float Edep_T6S = (intADC_T6S * calib_T6S);
        float Edep_T7S = (intADC_T7S * calib_T7S);
        float Edep_T8S = (intADC_T8S * calib_T8S);
        float Edep_T9S = (intADC_T9S * calib_T9S);

        float Edep_total_C =    Edep_T1C + Edep_T2C + Edep_T3C +
                                Edep_T4C + Edep_T5C + Edep_T6C +
                                Edep_T7C + Edep_T8C + Edep_T9C;

        float Edep_total_S =    Edep_T1S + Edep_T2S + Edep_T3S +
                                Edep_T4S + Edep_T5S + Edep_T6S +
                                Edep_T7S + Edep_T8S + Edep_T9S;

        float Edep_total_Comb = Edep_total_C + Edep_total_S;

        float scaledEdep_total_C = Edep_total_C * scaleFactor_C;
        float scaledEdep_total_S = Edep_total_S * scaleFactor_S;
        float scaledEdep_total_Comb = scaledEdep_total_C + scaledEdep_total_S;

        //////////////////////////////////////////////////////////////////////
        // Filling histograms before event selection
        //////////////////////////////////////////////////////////////////////
        hist_T1C_intADC->Fill(intADC_T1C);
        hist_T2C_intADC->Fill(intADC_T2C);
        hist_T3C_intADC->Fill(intADC_T3C);
        hist_T4C_intADC->Fill(intADC_T4C);
        hist_T5C_intADC->Fill(intADC_T5C);
        hist_T6C_intADC->Fill(intADC_T6C);
        hist_T7C_intADC->Fill(intADC_T7C);
        hist_T8C_intADC->Fill(intADC_T8C);
        hist_T9C_intADC->Fill(intADC_T9C);

        hist_T1S_intADC->Fill(intADC_T1S);
        hist_T2S_intADC->Fill(intADC_T2S);
        hist_T3S_intADC->Fill(intADC_T3S);
        hist_T4S_intADC->Fill(intADC_T4S);
        hist_T5S_intADC->Fill(intADC_T5S);
        hist_T6S_intADC->Fill(intADC_T6S);
        hist_T7S_intADC->Fill(intADC_T7S);
        hist_T8S_intADC->Fill(intADC_T8S);
        hist_T9S_intADC->Fill(intADC_T9S);

        hist_T1C_peakADC->Fill(peakADC_T1C);
        hist_T2C_peakADC->Fill(peakADC_T2C);
        hist_T3C_peakADC->Fill(peakADC_T3C);
        hist_T4C_peakADC->Fill(peakADC_T4C);
        hist_T5C_peakADC->Fill(peakADC_T5C);
        hist_T6C_peakADC->Fill(peakADC_T6C);
        hist_T7C_peakADC->Fill(peakADC_T7C);
        hist_T8C_peakADC->Fill(peakADC_T8C);
        hist_T9C_peakADC->Fill(peakADC_T9C);

        hist_T1S_peakADC->Fill(peakADC_T1S);
        hist_T2S_peakADC->Fill(peakADC_T2S);
        hist_T3S_peakADC->Fill(peakADC_T3S);
        hist_T4S_peakADC->Fill(peakADC_T4S);
        hist_T5S_peakADC->Fill(peakADC_T5S);
        hist_T6S_peakADC->Fill(peakADC_T6S);
        hist_T7S_peakADC->Fill(peakADC_T7S);
        hist_T8S_peakADC->Fill(peakADC_T8S);
        hist_T9S_peakADC->Fill(peakADC_T9S);

        hist_Edep_T1C->Fill(Edep_T1C);
        hist_Edep_T2C->Fill(Edep_T2C);
        hist_Edep_T3C->Fill(Edep_T3C);
        hist_Edep_T4C->Fill(Edep_T4C);
        hist_Edep_T5C->Fill(Edep_T5C);
        hist_Edep_T6C->Fill(Edep_T6C);
        hist_Edep_T7C->Fill(Edep_T7C);
        hist_Edep_T8C->Fill(Edep_T8C);
        hist_Edep_T9C->Fill(Edep_T9C);

        hist_Edep_T1S->Fill(Edep_T1S);
        hist_Edep_T2S->Fill(Edep_T2S);
        hist_Edep_T3S->Fill(Edep_T3S);
        hist_Edep_T4S->Fill(Edep_T4S);
        hist_Edep_T5S->Fill(Edep_T5S);
        hist_Edep_T6S->Fill(Edep_T6S);
        hist_Edep_T7S->Fill(Edep_T7S);
        hist_Edep_T8S->Fill(Edep_T8S);
        hist_Edep_T9S->Fill(Edep_T9S);

        hist_totalEdep_C->Fill(Edep_total_C);
        hist_totalEdep_S->Fill(Edep_total_S);
        hist_totalEdep_Comb->Fill(Edep_total_Comb);

        hist_totalEdepScaled_C->Fill(scaledEdep_total_C);
        hist_totalEdepScaled_S->Fill(scaledEdep_total_S);
        hist_totalEdepScaled_Comb->Fill(scaledEdep_total_Comb);

        hist_pos_WCX->Fill(pos_X);
        hist_pos_WCY->Fill(pos_Y);
        hist_WC->Fill(pos_X, pos_Y);
        hist_WC_fine->Fill(pos_X, pos_Y);

        //////////////////////////////////////////////////////////////////////
        // Event selection (Wire chamber only; hodoscope not available)
        //////////////////////////////////////////////////////////////////////
        bool passWCXcut = (pos_X > cut_WC_X_min && pos_X < cut_WC_X_max); // Require -5 mm < X position < 5 mm
        bool passWCYcut = (pos_Y > cut_WC_Y_min && pos_Y < cut_WC_Y_max); // Require -5 mm < Y position < 5 mm

        if (! (passWCXcut && passWCYcut) ) continue; // Fail event selection if not in the window

        //////////////////////////////////////////////////////////////////////
        // Filling histograms after event selection
        //////////////////////////////////////////////////////////////////////
        hist_T1C_intADC_after->Fill(intADC_T1C);
        hist_T2C_intADC_after->Fill(intADC_T2C);
        hist_T3C_intADC_after->Fill(intADC_T3C);
        hist_T4C_intADC_after->Fill(intADC_T4C);
        hist_T5C_intADC_after->Fill(intADC_T5C);
        hist_T6C_intADC_after->Fill(intADC_T6C);
        hist_T7C_intADC_after->Fill(intADC_T7C);
        hist_T8C_intADC_after->Fill(intADC_T8C);
        hist_T9C_intADC_after->Fill(intADC_T9C);

        hist_T1S_intADC_after->Fill(intADC_T1S);
        hist_T2S_intADC_after->Fill(intADC_T2S);
        hist_T3S_intADC_after->Fill(intADC_T3S);
        hist_T4S_intADC_after->Fill(intADC_T4S);
        hist_T5S_intADC_after->Fill(intADC_T5S);
        hist_T6S_intADC_after->Fill(intADC_T6S);
        hist_T7S_intADC_after->Fill(intADC_T7S);
        hist_T8S_intADC_after->Fill(intADC_T8S);
        hist_T9S_intADC_after->Fill(intADC_T9S);

        hist_T1C_peakADC_after->Fill(peakADC_T1C);
        hist_T2C_peakADC_after->Fill(peakADC_T2C);
        hist_T3C_peakADC_after->Fill(peakADC_T3C);
        hist_T4C_peakADC_after->Fill(peakADC_T4C);
        hist_T5C_peakADC_after->Fill(peakADC_T5C);
        hist_T6C_peakADC_after->Fill(peakADC_T6C);
        hist_T7C_peakADC_after->Fill(peakADC_T7C);
        hist_T8C_peakADC_after->Fill(peakADC_T8C);
        hist_T9C_peakADC_after->Fill(peakADC_T9C);

        hist_T1S_peakADC_after->Fill(peakADC_T1S);
        hist_T2S_peakADC_after->Fill(peakADC_T2S);
        hist_T3S_peakADC_after->Fill(peakADC_T3S);
        hist_T4S_peakADC_after->Fill(peakADC_T4S);
        hist_T5S_peakADC_after->Fill(peakADC_T5S);
        hist_T6S_peakADC_after->Fill(peakADC_T6S);
        hist_T7S_peakADC_after->Fill(peakADC_T7S);
        hist_T8S_peakADC_after->Fill(peakADC_T8S);
        hist_T9S_peakADC_after->Fill(peakADC_T9S);

        hist_Edep_T1C_after->Fill(Edep_T1C);
        hist_Edep_T2C_after->Fill(Edep_T2C);
        hist_Edep_T3C_after->Fill(Edep_T3C);
        hist_Edep_T4C_after->Fill(Edep_T4C);
        hist_Edep_T5C_after->Fill(Edep_T5C);
        hist_Edep_T6C_after->Fill(Edep_T6C);
        hist_Edep_T7C_after->Fill(Edep_T7C);
        hist_Edep_T8C_after->Fill(Edep_T8C);
        hist_Edep_T9C_after->Fill(Edep_T9C);

        hist_Edep_T1S_after->Fill(Edep_T1S);
        hist_Edep_T2S_after->Fill(Edep_T2S);
        hist_Edep_T3S_after->Fill(Edep_T3S);
        hist_Edep_T4S_after->Fill(Edep_T4S);
        hist_Edep_T5S_after->Fill(Edep_T5S);
        hist_Edep_T6S_after->Fill(Edep_T6S);
        hist_Edep_T7S_after->Fill(Edep_T7S);
        hist_Edep_T8S_after->Fill(Edep_T8S);
        hist_Edep_T9S_after->Fill(Edep_T9S);

        hist_totalEdep_C_after->Fill(Edep_total_C);
        hist_totalEdep_S_after->Fill(Edep_total_S);
        hist_totalEdep_Comb_after->Fill(Edep_total_Comb);

        hist_totalEdepScaled_C_after->Fill(scaledEdep_total_C);
        hist_totalEdepScaled_S_after->Fill(scaledEdep_total_S);
        hist_totalEdepScaled_Comb_after->Fill(scaledEdep_total_Comb);

        hist_pos_WCX_after->Fill(pos_X);
        hist_pos_WCY_after->Fill(pos_Y);
        hist_WC_after->Fill(pos_X, pos_Y);
        hist_WC_fine_after->Fill(pos_X, pos_Y);
    }

    //////////////////////////////////////////////////////////////////////
    // Fitting functions
    //////////////////////////////////////////////////////////////////////
    // Fit function for calculating Scale factor
    TF1 *fit_scale_C = new TF1("fit_scale_C","gaus",0,10);
    fit_scale_C->SetParameters(1,0,1);
    hist_totalEdep_C_after->Fit("fit_scale_C");

    TF1 *fit_scale_S = new TF1("fit_scale_S","gaus",0,10);
    fit_scale_S->SetParameters(1,0,1);
    hist_totalEdep_S_after->Fit("fit_scale_S");

    // Fit function for energy resolution fit
    TF1 *fit_reso_C = new TF1("fit_reso_C","gaus",0,10);
    fit_reso_C->SetParameters(1,0,1);
    hist_totalEdepScaled_C_after->Fit("fit_reso_C");

    TF1 *fit_reso_S = new TF1("fit_reso_S","gaus",0,10);
    fit_reso_S->SetParameters(1,0,1);
    hist_totalEdepScaled_S_after->Fit("fit_reso_S");

    TF1 *fit_reso_Comb = new TF1("fit_reso_Comb","gaus",0,20);
    fit_reso_Comb->SetParameters(1,0,1);
    hist_totalEdepScaled_Comb_after->Fit("fit_reso_Comb");

    //////////////////////////////////////////////////////////////////////
    // Output file
    //////////////////////////////////////////////////////////////////////
    std::string outFile = "./Calib/Calib_Run_" + std::to_string(fRunNum) + ".root";
    TFile* outputRoot = new TFile(outFile.c_str(), "RECREATE");
    outputRoot->cd();

    //////////////////////////////////////////////////////////////////////
    // Writing fit functions
    //////////////////////////////////////////////////////////////////////
    fit_scale_C->Write();
    fit_scale_S->Write();

    fit_reso_C->Write();
    fit_reso_S->Write();
    fit_reso_Comb->Write();

    //////////////////////////////////////////////////////////////////////
    // Writing histograms before event selection
    //////////////////////////////////////////////////////////////////////
    hist_T1C_intADC->Write();
    hist_T2C_intADC->Write();
    hist_T3C_intADC->Write();
    hist_T4C_intADC->Write();
    hist_T5C_intADC->Write();
    hist_T6C_intADC->Write();
    hist_T7C_intADC->Write();
    hist_T8C_intADC->Write();
    hist_T9C_intADC->Write();

    hist_T1S_intADC->Write();
    hist_T2S_intADC->Write();
    hist_T3S_intADC->Write();
    hist_T4S_intADC->Write();
    hist_T5S_intADC->Write();
    hist_T6S_intADC->Write();
    hist_T7S_intADC->Write();
    hist_T8S_intADC->Write();
    hist_T9S_intADC->Write();

    hist_T1C_peakADC->Write();
    hist_T2C_peakADC->Write();
    hist_T3C_peakADC->Write();
    hist_T4C_peakADC->Write();
    hist_T5C_peakADC->Write();
    hist_T6C_peakADC->Write();
    hist_T7C_peakADC->Write();
    hist_T8C_peakADC->Write();
    hist_T9C_peakADC->Write();

    hist_T1S_peakADC->Write();
    hist_T2S_peakADC->Write();
    hist_T3S_peakADC->Write();
    hist_T4S_peakADC->Write();
    hist_T5S_peakADC->Write();
    hist_T6S_peakADC->Write();
    hist_T7S_peakADC->Write();
    hist_T8S_peakADC->Write();
    hist_T9S_peakADC->Write();

    hist_Edep_T1C->Write();
    hist_Edep_T2C->Write();
    hist_Edep_T3C->Write();
    hist_Edep_T4C->Write();
    hist_Edep_T5C->Write();
    hist_Edep_T6C->Write();
    hist_Edep_T7C->Write();
    hist_Edep_T8C->Write();
    hist_Edep_T9C->Write();

    hist_Edep_T1S->Write();
    hist_Edep_T2S->Write();
    hist_Edep_T3S->Write();
    hist_Edep_T4S->Write();
    hist_Edep_T5S->Write();
    hist_Edep_T6S->Write();
    hist_Edep_T7S->Write();
    hist_Edep_T8S->Write();
    hist_Edep_T9S->Write();

    hist_totalEdep_C->Write();
    hist_totalEdep_S->Write();
    hist_totalEdep_Comb->Write();

    hist_totalEdepScaled_C->Write();
    hist_totalEdepScaled_S->Write();
    hist_totalEdepScaled_Comb->Write();

    hist_pos_WCX->Write();
    hist_pos_WCY->Write();
    hist_WC->Write();
    hist_WC_fine->Write();

    //////////////////////////////////////////////////////////////////////
    // Writing histograms after event selection
    //////////////////////////////////////////////////////////////////////
    hist_T1C_intADC_after->Write();
    hist_T2C_intADC_after->Write();
    hist_T3C_intADC_after->Write();
    hist_T4C_intADC_after->Write();
    hist_T5C_intADC_after->Write();
    hist_T6C_intADC_after->Write();
    hist_T7C_intADC_after->Write();
    hist_T8C_intADC_after->Write();
    hist_T9C_intADC_after->Write();

    hist_T1S_intADC_after->Write();
    hist_T2S_intADC_after->Write();
    hist_T3S_intADC_after->Write();
    hist_T4S_intADC_after->Write();
    hist_T5S_intADC_after->Write();
    hist_T6S_intADC_after->Write();
    hist_T7S_intADC_after->Write();
    hist_T8S_intADC_after->Write();
    hist_T9S_intADC_after->Write();

    hist_T1C_peakADC_after->Write();
    hist_T2C_peakADC_after->Write();
    hist_T3C_peakADC_after->Write();
    hist_T4C_peakADC_after->Write();
    hist_T5C_peakADC_after->Write();
    hist_T6C_peakADC_after->Write();
    hist_T7C_peakADC_after->Write();
    hist_T8C_peakADC_after->Write();
    hist_T9C_peakADC_after->Write();

    hist_T1S_peakADC_after->Write();
    hist_T2S_peakADC_after->Write();
    hist_T3S_peakADC_after->Write();
    hist_T4S_peakADC_after->Write();
    hist_T5S_peakADC_after->Write();
    hist_T6S_peakADC_after->Write();
    hist_T7S_peakADC_after->Write();
    hist_T8S_peakADC_after->Write();
    hist_T9S_peakADC_after->Write();

    hist_Edep_T1C_after->Write();
    hist_Edep_T2C_after->Write();
    hist_Edep_T3C_after->Write();
    hist_Edep_T4C_after->Write();
    hist_Edep_T5C_after->Write();
    hist_Edep_T6C_after->Write();
    hist_Edep_T7C_after->Write();
    hist_Edep_T8C_after->Write();
    hist_Edep_T9C_after->Write();

    hist_Edep_T1S_after->Write();
    hist_Edep_T2S_after->Write();
    hist_Edep_T3S_after->Write();
    hist_Edep_T4S_after->Write();
    hist_Edep_T5S_after->Write();
    hist_Edep_T6S_after->Write();
    hist_Edep_T7S_after->Write();
    hist_Edep_T8S_after->Write();
    hist_Edep_T9S_after->Write();

    hist_totalEdep_C_after->Write();
    hist_totalEdep_S_after->Write();
    hist_totalEdep_Comb_after->Write();

    hist_totalEdepScaled_C_after->Write();
    hist_totalEdepScaled_S_after->Write();
    hist_totalEdepScaled_Comb_after->Write();

    hist_pos_WCX_after->Write();
    hist_pos_WCY_after->Write();
    hist_WC_after->Write();
    hist_WC_fine_after->Write();

    outputRoot->Close();
}
