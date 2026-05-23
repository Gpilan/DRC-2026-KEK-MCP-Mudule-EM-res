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
    util.LoadMapping("../../mapping/KEK_TB2026May_PMT.root");
    // util.LoadMapping("../../mapping/KEK_TB2026May_MCP.root");

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

    TBcid cid_trg1 = util.GetCID("Trg1");
    TBcid cid_trg2 = util.GetCID("Trg2");

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

    TH1F* hist_trg1 = new TH1F("Trigger 1", ";bin;ADC", 1000, 0, 1000);
    TH1F* hist_trg2 = new TH1F("Trigger 2", ";bin;ADC", 1000, 0, 1000);

    // Change {8, 9, 13} to actual set of MIDs
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

    hist_trg1->Write();
    hist_trg2->Write();

    outputRoot->Close();
}
