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

    fs::path dir("./Avg");   
    if (!(fs::exists(dir))) fs::create_directory(dir);

    // initialize the utility class
    TButility util = TButility();
    util.LoadMapping("../mapping/mapping_KEK.root");
    
    // prepare CIDs that we want to use (CID = Channel ID)
    // 3x3 tower CIDs
    TBcid cid_T1_C = util.GetCID("T1-C");
    TBcid cid_T2_C = util.GetCID("T2-C");
    TBcid cid_T3_C = util.GetCID("T3-C");
    TBcid cid_T4_C = util.GetCID("T4-C");
    TBcid cid_T5_C = util.GetCID("T5-C");
    TBcid cid_T6_C = util.GetCID("T6-C");
    TBcid cid_T7_C = util.GetCID("T7-C");
    TBcid cid_T8_C = util.GetCID("T8-C");
    TBcid cid_T9_C = util.GetCID("T9-C");

    TBcid cid_T1_S = util.GetCID("T1-S");
    TBcid cid_T2_S = util.GetCID("T2-S");
    TBcid cid_T3_S = util.GetCID("T3-S");
    TBcid cid_T4_S = util.GetCID("T4-S");
    TBcid cid_T5_S = util.GetCID("T5-S");
    TBcid cid_T6_S = util.GetCID("T6-S");
    TBcid cid_T7_S = util.GetCID("T7-S");
    TBcid cid_T8_S = util.GetCID("T8-S");
    TBcid cid_T9_S = util.GetCID("T9-S");
    // Aux. detectors
    // Wire chamber
    TBcid cid_WCX = util.GetCID("WCX"); // Wire chamber X
    TBcid cid_WCY = util.GetCID("WCY"); // Wire chamber Y

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

    TBcid cid_trg1 = util.GetCID("trg1");
    TBcid cid_trg2 = util.GetCID("trg2");

    // prepare the histograms wa want to draw
    TH1F* hist_T1_C = new TH1F("T1_C" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T2_C = new TH1F("T2_C" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T3_C = new TH1F("T3_C" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T4_C = new TH1F("T4_C" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T5_C = new TH1F("T5_C" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T6_C = new TH1F("T6_C" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T7_C = new TH1F("T7_C" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T8_C = new TH1F("T8_C" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T9_C = new TH1F("T9_C" , ";bin;nEvents", 1000, 0, 1000);

    TH1F* hist_T1_S = new TH1F("T1_S" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T2_S = new TH1F("T2_S" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T3_S = new TH1F("T3_S" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T4_S = new TH1F("T4_S" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T5_S = new TH1F("T5_S" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T6_S = new TH1F("T6_S" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T7_S = new TH1F("T7_S" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T8_S = new TH1F("T8_S" , ";bin;nEvents", 1000, 0, 1000);
    TH1F* hist_T9_S = new TH1F("T9_S" , ";bin;nEvents", 1000, 0, 1000);

    TH1F* hist_WCX = new TH1F("Wire chamber X", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_WCY = new TH1F("Wire chamber Y", ";bin;ADC", 1000, 0, 1000);

    TH1F* hist_HX1 = new TH1F("Hodoscope X1", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX2 = new TH1F("Hodoscope X2", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX3 = new TH1F("Hodoscope X3", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX4 = new TH1F("Hodoscope X4", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX5 = new TH1F("Hodoscope X5", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX6 = new TH1F("Hodoscope X6", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX7 = new TH1F("Hodoscope X7", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX8 = new TH1F("Hodoscope X8", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX9 = new TH1F("Hodoscope X9", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX10 = new TH1F("Hodoscope X10", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX11 = new TH1F("Hodoscope X11", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX12 = new TH1F("Hodoscope X12", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX13 = new TH1F("Hodoscope X13", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX14 = new TH1F("Hodoscope X14", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX15 = new TH1F("Hodoscope X15", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HX16 = new TH1F("Hodoscope X16", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY1 = new TH1F("Hodoscope Y1", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY2 = new TH1F("Hodoscope Y2", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY3 = new TH1F("Hodoscope Y3", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY4 = new TH1F("Hodoscope Y4", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY5 = new TH1F("Hodoscope Y5", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY6 = new TH1F("Hodoscope Y6", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY7 = new TH1F("Hodoscope Y7", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY8 = new TH1F("Hodoscope Y8", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY9 = new TH1F("Hodoscope Y9", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY10 = new TH1F("Hodoscope Y10", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY11 = new TH1F("Hodoscope Y11", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY12 = new TH1F("Hodoscope Y12", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY13 = new TH1F("Hodoscope Y13", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY14 = new TH1F("Hodoscope Y14", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY15 = new TH1F("Hodoscope Y15", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_HY16 = new TH1F("Hodoscope Y16", ";bin;ADC", 1000, 0, 1000);

    TH1F* hist_trg1 = new TH1F("Trigger 1", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_trg2 = new TH1F("Trigger 2", ";bin;ADC", 1000, 0, 1000);

    // Change {8, 9, 13} to actual set of MIDs
    TBread<TBwaveform> readerWave = TBread<TBwaveform>(fRunNum, fMaxEvent, -1, false, "/pnfs/knu.ac.kr/data/cms/store/user/sungwon/2026_May_KEK_DRC_TB_Data", {8, 9, 13});

    // Set Maximum event
    if (fMaxEvent == -1 || fMaxEvent > readerWave.GetMaxEvent())
        fMaxEvent = readerWave.GetMaxEvent();

    // Evt Loop
    for (int iEvt = 0; iEvt < fMaxEvent; iEvt++) {
        printProgress(iEvt, fMaxEvent);
        TBevt<TBwaveform> anEvt = readerWave.GetAnEvent();

        std::vector<short> wave_T1_C = anEvt.GetData(cid_T1_C).waveform();
        std::vector<short> wave_T2_C = anEvt.GetData(cid_T2_C).waveform();
        std::vector<short> wave_T3_C = anEvt.GetData(cid_T3_C).waveform();
        std::vector<short> wave_T4_C = anEvt.GetData(cid_T4_C).waveform();
        std::vector<short> wave_T5_C = anEvt.GetData(cid_T5_C).waveform();
        std::vector<short> wave_T6_C = anEvt.GetData(cid_T6_C).waveform();
        std::vector<short> wave_T7_C = anEvt.GetData(cid_T7_C).waveform();
        std::vector<short> wave_T8_C = anEvt.GetData(cid_T8_C).waveform();
        std::vector<short> wave_T9_C = anEvt.GetData(cid_T9_C).waveform();

        std::vector<short> wave_T1_S = anEvt.GetData(cid_T1_S).waveform();
        std::vector<short> wave_T2_S = anEvt.GetData(cid_T2_S).waveform();
        std::vector<short> wave_T3_S = anEvt.GetData(cid_T3_S).waveform();
        std::vector<short> wave_T4_S = anEvt.GetData(cid_T4_S).waveform();
        std::vector<short> wave_T5_S = anEvt.GetData(cid_T5_S).waveform();
        std::vector<short> wave_T6_S = anEvt.GetData(cid_T6_S).waveform();
        std::vector<short> wave_T7_S = anEvt.GetData(cid_T7_S).waveform();
        std::vector<short> wave_T8_S = anEvt.GetData(cid_T8_S).waveform();
        std::vector<short> wave_T9_S = anEvt.GetData(cid_T9_S).waveform();

        // For trigger
        std::vector<short> wave_trg1 = (anEvt.GetData(cid_trg1)).waveform();
        std::vector<short> wave_trg2 = (anEvt.GetData(cid_trg2)).waveform();

        // Get waveform for DWCs
        std::vector<short> wave_WCX = (anEvt.GetData(cid_WCX)).waveform();
        std::vector<short> wave_WCY = (anEvt.GetData(cid_WCY)).waveform();

        std::vector<short> wave_HX1 = (anEvt.GetData(cid_HX1)).waveform();
        std::vector<short> wave_HX2 = (anEvt.GetData(cid_HX2)).waveform();
        std::vector<short> wave_HX3 = (anEvt.GetData(cid_HX3)).waveform();
        std::vector<short> wave_HX4 = (anEvt.GetData(cid_HX4)).waveform();
        std::vector<short> wave_HX5 = (anEvt.GetData(cid_HX5)).waveform();
        std::vector<short> wave_HX6 = (anEvt.GetData(cid_HX6)).waveform();
        std::vector<short> wave_HX7 = (anEvt.GetData(cid_HX7)).waveform();
        std::vector<short> wave_HX8 = (anEvt.GetData(cid_HX8)).waveform();
        std::vector<short> wave_HX9 = (anEvt.GetData(cid_HX9)).waveform();
        std::vector<short> wave_HX10 = (anEvt.GetData(cid_HX10)).waveform();
        std::vector<short> wave_HX11 = (anEvt.GetData(cid_HX11)).waveform();
        std::vector<short> wave_HX12 = (anEvt.GetData(cid_HX12)).waveform();
        std::vector<short> wave_HX13 = (anEvt.GetData(cid_HX13)).waveform();
        std::vector<short> wave_HX14 = (anEvt.GetData(cid_HX14)).waveform();
        std::vector<short> wave_HX15 = (anEvt.GetData(cid_HX15)).waveform();    
        std::vector<short> wave_HX16 = (anEvt.GetData(cid_HX16)).waveform();

        std::vector<short> wave_HY1 = (anEvt.GetData(cid_HY1)).waveform();
        std::vector<short> wave_HY2 = (anEvt.GetData(cid_HY2)).waveform();
        std::vector<short> wave_HY3 = (anEvt.GetData(cid_HY3)).waveform();
        std::vector<short> wave_HY4 = (anEvt.GetData(cid_HY4)).waveform();
        std::vector<short> wave_HY5 = (anEvt.GetData(cid_HY5)).waveform();
        std::vector<short> wave_HY6 = (anEvt.GetData(cid_HY6)).waveform();
        std::vector<short> wave_HY7 = (anEvt.GetData(cid_HY7)).waveform();
        std::vector<short> wave_HY8 = (anEvt.GetData(cid_HY8)).waveform();
        std::vector<short> wave_HY9 = (anEvt.GetData(cid_HY9)).waveform();
        std::vector<short> wave_HY10 = (anEvt.GetData(cid_HY10)).waveform();
        std::vector<short> wave_HY11 = (anEvt.GetData(cid_HY11)).waveform();
        std::vector<short> wave_HY12 = (anEvt.GetData(cid_HY12)).waveform();
        std::vector<short> wave_HY13 = (anEvt.GetData(cid_HY13)).waveform();
        std::vector<short> wave_HY14 = (anEvt.GetData(cid_HY14)).waveform();
        std::vector<short> wave_HY15 = (anEvt.GetData(cid_HY15)).waveform();
        std::vector<short> wave_HY16 = (anEvt.GetData(cid_HY16)).waveform();

        for (int bin = 1; bin < 1001; bin++) {      

            hist_T1_C->Fill(bin, (float) (wave_T1_C.at(bin)) / (float) (fMaxEvent) );
            hist_T2_C->Fill(bin, (float) (wave_T2_C.at(bin)) / (float) (fMaxEvent) );
            hist_T3_C->Fill(bin, (float) (wave_T3_C.at(bin)) / (float) (fMaxEvent) );
            hist_T4_C->Fill(bin, (float) (wave_T4_C.at(bin)) / (float) (fMaxEvent) );
            hist_T5_C->Fill(bin, (float) (wave_T5_C.at(bin)) / (float) (fMaxEvent) );
            hist_T6_C->Fill(bin, (float) (wave_T6_C.at(bin)) / (float) (fMaxEvent) );
            hist_T7_C->Fill(bin, (float) (wave_T7_C.at(bin)) / (float) (fMaxEvent) );
            hist_T8_C->Fill(bin, (float) (wave_T8_C.at(bin)) / (float) (fMaxEvent) );
            hist_T9_C->Fill(bin, (float) (wave_T9_C.at(bin)) / (float) (fMaxEvent) );

            hist_T1_S->Fill(bin, (float) (wave_T1_S.at(bin)) / (float) (fMaxEvent) );
            hist_T2_S->Fill(bin, (float) (wave_T2_S.at(bin)) / (float) (fMaxEvent) );
            hist_T3_S->Fill(bin, (float) (wave_T3_S.at(bin)) / (float) (fMaxEvent) );
            hist_T4_S->Fill(bin, (float) (wave_T4_S.at(bin)) / (float) (fMaxEvent) );
            hist_T5_S->Fill(bin, (float) (wave_T5_S.at(bin)) / (float) (fMaxEvent) );
            hist_T6_S->Fill(bin, (float) (wave_T6_S.at(bin)) / (float) (fMaxEvent) );
            hist_T7_S->Fill(bin, (float) (wave_T7_S.at(bin)) / (float) (fMaxEvent) );
            hist_T8_S->Fill(bin, (float) (wave_T8_S.at(bin)) / (float) (fMaxEvent) );
            hist_T9_S->Fill(bin, (float) (wave_T9_S.at(bin)) / (float) (fMaxEvent) );

            hist_WCX->Fill(bin, (float) (wave_WCX.at(bin)) / (float) (fMaxEvent) );
            hist_WCY->Fill(bin, (float) (wave_WCY.at(bin)) / (float) (fMaxEvent) );

            hist_HX1->Fill(bin, (float) (wave_HX1.at(bin)) / (float) (fMaxEvent) );
            hist_HX2->Fill(bin, (float) (wave_HX2.at(bin)) / (float) (fMaxEvent) );
            hist_HX3->Fill(bin, (float) (wave_HX3.at(bin)) / (float) (fMaxEvent) );
            hist_HX4->Fill(bin, (float) (wave_HX4.at(bin)) / (float) (fMaxEvent) );
            hist_HX5->Fill(bin, (float) (wave_HX5.at(bin)) / (float) (fMaxEvent) );
            hist_HX6->Fill(bin, (float) (wave_HX6.at(bin)) / (float) (fMaxEvent) );
            hist_HX7->Fill(bin, (float) (wave_HX7.at(bin)) / (float) (fMaxEvent) );
            hist_HX8->Fill(bin, (float) (wave_HX8.at(bin)) / (float) (fMaxEvent) );
            hist_HX9->Fill(bin, (float) (wave_HX9.at(bin)) / (float) (fMaxEvent) );
            hist_HX10->Fill(bin, (float) (wave_HX10.at(bin)) / (float) (fMaxEvent) );
            hist_HX11->Fill(bin, (float) (wave_HX11.at(bin)) / (float) (fMaxEvent) );
            hist_HX12->Fill(bin, (float) (wave_HX12.at(bin)) / (float) (fMaxEvent) );
            hist_HX13->Fill(bin, (float) (wave_HX13.at(bin)) / (float) (fMaxEvent) );
            hist_HX14->Fill(bin, (float) (wave_HX14.at(bin)) / (float) (fMaxEvent) );
            hist_HX15->Fill(bin, (float) (wave_HX15.at(bin)) / (float) (fMaxEvent) );
            hist_HX16->Fill(bin, (float) (wave_HX16.at(bin)) / (float) (fMaxEvent) );

            hist_HY1->Fill(bin, (float) (wave_HY1.at(bin)) / (float) (fMaxEvent) );
            hist_HY2->Fill(bin, (float) (wave_HY2.at(bin)) / (float) (fMaxEvent) );
            hist_HY3->Fill(bin, (float) (wave_HY3.at(bin)) / (float) (fMaxEvent) );
            hist_HY4->Fill(bin, (float) (wave_HY4.at(bin)) / (float) (fMaxEvent) );
            hist_HY5->Fill(bin, (float) (wave_HY5.at(bin)) / (float) (fMaxEvent) );
            hist_HY6->Fill(bin, (float) (wave_HY6.at(bin)) / (float) (fMaxEvent) );
            hist_HY7->Fill(bin, (float) (wave_HY7.at(bin)) / (float) (fMaxEvent) );
            hist_HY8->Fill(bin, (float) (wave_HY8.at(bin)) / (float) (fMaxEvent) );
            hist_HY9->Fill(bin, (float) (wave_HY9.at(bin)) / (float) (fMaxEvent) );
            hist_HY10->Fill(bin, (float) (wave_HY10.at(bin)) / (float) (fMaxEvent) );
            hist_HY11->Fill(bin, (float) (wave_HY11.at(bin)) / (float) (fMaxEvent) );
            hist_HY12->Fill(bin, (float) (wave_HY12.at(bin)) / (float) (fMaxEvent) );
            hist_HY13->Fill(bin, (float) (wave_HY13.at(bin)) / (float) (fMaxEvent) );
            hist_HY14->Fill(bin, (float) (wave_HY14.at(bin)) / (float) (fMaxEvent) );
            hist_HY15->Fill(bin, (float) (wave_HY15.at(bin)) / (float) (fMaxEvent) );
            hist_HY16->Fill(bin, (float) (wave_HY16.at(bin)) / (float) (fMaxEvent) );

            hist_trg1->Fill(bin, (float) (wave_trg1.at(bin)) / (float) (fMaxEvent) );
            hist_trg2->Fill(bin, (float) (wave_trg2.at(bin)) / (float) (fMaxEvent) );
        }
    }

    std::string outFile = "./Avg/Avg_Run_" + std::to_string(fRunNum) + ".root";
    TFile* outputRoot = new TFile(outFile.c_str(), "RECREATE");
    outputRoot->cd();

    hist_T1_C->Write();
    hist_T2_C->Write();
    hist_T3_C->Write();
    hist_T4_C->Write();
    hist_T5_C->Write();
    hist_T6_C->Write();
    hist_T7_C->Write();
    hist_T8_C->Write();
    hist_T9_C->Write();

    hist_T1_S->Write();
    hist_T2_S->Write();
    hist_T3_S->Write();
    hist_T4_S->Write();
    hist_T5_S->Write();
    hist_T6_S->Write();
    hist_T7_S->Write();
    hist_T8_S->Write();
    hist_T9_S->Write();

    hist_WCX->Write();
    hist_WCY->Write();  

    hist_HX1->Write();
    hist_HX2->Write();
    hist_HX3->Write();
    hist_HX4->Write();
    hist_HX5->Write();
    hist_HX6->Write();
    hist_HX7->Write();
    hist_HX8->Write();
    hist_HX9->Write();
    hist_HX10->Write();
    hist_HX11->Write();
    hist_HX12->Write();
    hist_HX13->Write();
    hist_HX14->Write();
    hist_HX15->Write();
    hist_HX16->Write();

    hist_HY1->Write();
    hist_HY2->Write();
    hist_HY3->Write();
    hist_HY4->Write();
    hist_HY5->Write();
    hist_HY6->Write();
    hist_HY7->Write();
    hist_HY8->Write();
    hist_HY9->Write();
    hist_HY10->Write();
    hist_HY11->Write();
    hist_HY12->Write();
    hist_HY13->Write();
    hist_HY14->Write();
    hist_HY15->Write();
    hist_HY16->Write();

    hist_trg1->Write();
    hist_trg2->Write();

    outputRoot->Close();
}
