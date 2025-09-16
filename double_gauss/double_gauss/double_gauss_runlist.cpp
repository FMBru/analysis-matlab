#include <TColor.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <sstream>


using namespace std;


double fitDoubleGaussian(double *x, double *par) {

    double mean, sigmaCore, sigmaTail, a, scale;
    double coreGauss, tailGauss;

    mean = par[0];
    sigmaCore = par[1];
    sigmaTail = par[2];
    a = par[3];
    scale = par[4];

    coreGauss = a * TMath::Exp(- TMath::Power(x[0] - mean, 2) / (2 * TMath::Power(sigmaCore, 2)));
    tailGauss = (1 - a) * TMath::Exp(- TMath::Power(x[0] - mean, 2) / (2 * TMath::Power(sigmaTail, 2)));

    return (scale * (coreGauss + tailGauss));

}

double getSigma(double a, double sigmaCore, double sigmaTail) {

    return TMath::Sqrt(a * TMath::Power(sigmaCore, 2) + (1 - a) * TMath::Power(sigmaTail, 2));

}

double getSigmaError(double a, double aError, double sigmaCore,
                     double sigmaCoreError, double sigmaTail, double sigmaTailError) {

    double aDerivative = (TMath::Power(sigmaCore, 2) - TMath::Power(sigmaTail, 2)) / (2 * TMath::Sqrt(a * (TMath::Power(sigmaCore, 2) - TMath::Power(sigmaTail, 2)) + TMath::Power(sigmaTail, 2)));
    double sigmaCoreDerivative = (a * sigmaCore) / (TMath::Sqrt(a * (TMath::Power(sigmaCore, 2) - TMath::Power(sigmaTail, 2)) + TMath::Power(sigmaTail, 2)));
    double sigmaTailDerivative = (- (a - 1) * sigmaTail) / (TMath::Sqrt(a * (TMath::Power(sigmaCore, 2) - TMath::Power(sigmaTail, 2)) + TMath::Power(sigmaTail, 2)));

    double aPart = TMath::Power(aDerivative * aError, 2);
    double sigmaCorePart = TMath::Power(sigmaCoreDerivative * sigmaCoreError, 2);
    double sigmaTailPart = TMath::Power(sigmaTailDerivative * sigmaTailError, 2);

    return TMath::Sqrt(aPart + sigmaCorePart + sigmaTailPart);

}

double myGaus(double *x, double *par) {

    double scale, mean, sigma;

    mean = par[1];
    sigma = par[2];
    scale = par[0];

    return scale * TMath::Exp(- TMath::Power(x[0] - mean, 2) / (2 * TMath::Power(sigma, 2)));

}


void double_gauss(const string& inputPath, const string& outputPath, const string& filename, int nbins) {
    
    // Construct full file paths
    string inputFile = inputPath + "/" + filename;
    string outputFileTxt = outputPath + "/" + filename.substr(0, 6) + "_ROOTfit.txt";
    string outputFilePdf = outputPath + "/" + filename.substr(0, filename.find_last_of('.')) + "_ROOTfit.pdf";
    
	double g1_min = 999999.0;
	double g1_max = -999999.0;
    std::vector<double> data_in;
	TCanvas *c1 = new TCanvas("c1","spectrum");
    //c1->SetLogy(1);
    c1->SetFillColor(0);
    c1->SetGrid(); 
    gStyle->SetOptStat(1111);
    //gStyle->SetOptFit(111);

	FILE *fp = fopen(inputFile.c_str(),"r");
	double data;
	int bin;
	int data_n=0;
	while (!feof(fp))
	{
        [[maybe_unused]] auto x = fscanf(fp,"%lf",&data);
		//h1->Fill(data);
		if (data>g1_max) g1_max = data;
		if (data<g1_min) g1_min = data;
		data_in.push_back(data);
		data_n++;
	}
	fclose(fp);
	

	TH1F *h1 = new TH1F("h1", "Time resolution",nbins, g1_min, g1_max);
	for(int i=0;i<data_n;i++)
	{
		h1->Fill(data_in.at(i));
	}
	h1->GetXaxis()->SetTitle("Time difference: PICOSEC vs reference, ns");
	h1->GetYaxis()->SetTitle("Entries");
	h1->Rebin(1);
	
    /* Fit using single gaussian*/
	TF1 *g1    = new TF1("g1","gaus",g1_min,g1_max);
	Double_t par[3];
	g1->SetLineWidth(2);
	g1->SetLineColor(kGreen);
	h1->Fit(g1,"L");
	g1->GetParameters(&par[0]);
	h1->Draw();

	TF1 *fit = new TF1("fit", fitDoubleGaussian, g1_min, g1_max, 5);
	fit->SetParameters(h1->GetMean(), h1->GetRMS() * 0.5, h1->GetRMS() * 1.5, 0.8, h1->GetBinContent(h1->GetMaximumBin()));
	fit->SetParNames("mean", "sigmaCore", "sigmaTail", "a", "scale");
	fit->SetLineColor(kRed);
    h1->Fit("fit", "L");

	int entries = h1->GetEntries();
	double stddev = h1->GetStdDev();
	double chi2 = fit->GetChisquare();
	double ndf = fit->GetNDF();
    double mean = fit->GetParameter("mean");
	double meanError = fit->GetParError(0);
	double sigmaCore = fit->GetParameter("sigmaCore");
	double sigmaCoreError = fit->GetParError(1);
	double sigmaTail = fit->GetParameter("sigmaTail");
	double sigmaTailError = fit->GetParError(2);
	double a = fit->GetParameter("a");
	double aError = fit->GetParError(3);
	double scale = fit->GetParameter("scale");
	double scaleError = fit->GetParError(4);
	
	double sigma = getSigma(a, sigmaCore, sigmaTail);
	double sigmaError = getSigmaError(a, aError, sigmaCore, sigmaCoreError, sigmaTail, sigmaTailError);
// 	cout << "SIGMA: "<< sigma*1000 << " +- " << sigmaError *1000<< endl;
	
// 	cout << stddev << endl;

	TF1 *fitCore = new TF1("fitCore", myGaus, g1_min, g1_max,3);
	fitCore->SetParameters(scale * a, mean, sigmaCore);
	fitCore->SetLineColor(kViolet);
	fitCore->SetLineStyle(kDashed);
	fitCore->Draw("same");

	TF1 *fitTail = new TF1("fitTail", myGaus, g1_min, g1_max,3);
	fitTail->SetParameters(scale * (1.0 - a), mean, sigmaTail);
	fitTail->SetLineColor(kMagenta);
	fitTail->SetLineStyle(kDashed);
	fitTail->Draw("same");


	g1->Draw("same");
	fit->Draw("same");

	// Build and Draw a legend
    TLegend leg(.1,.7,.4,.9,"Double Gauss");
    leg.SetFillColor(0);

//     TLatex ltx = TLatex();
//     ltx.SetTextSize(0.035);
//     char ltxText[256];
//     snprintf(ltxText, sizeof(ltxText), "#mu_{Single} = %.3f #pm %.3f ns", par[1], g1->GetParError(1));
//     leg.AddEntry(g1, ltxText);

    TLatex ltx1 = TLatex();
    ltx1.SetTextSize(0.035);
    char ltx1Text[256];
    snprintf(ltx1Text, sizeof(ltx1Text), "#sigma_{Single} = %.1f #pm %.1f ps", par[2]*1000, g1->GetParError(2)*1000);
    leg.AddEntry(g1, ltx1Text);

    TLatex ltx2 = TLatex();
    ltx2.SetTextSize(0.035);
    char ltx2Text[256];
    snprintf(ltx2Text, sizeof(ltx2Text), "#sigma_{Comb} = %.1f #pm %.1f ps", sigma*1000, sigmaError*1000);
    leg.AddEntry(fit, ltx2Text);

//     TLatex ltx3 = TLatex();
//     ltx3.SetTextSize(0.035);
//     char ltx3Text[256];
//     snprintf(ltx3Text, sizeof(ltx3Text), "#mu_{Comb} = %.3f #pm %.3f ns", mean, meanError);
//     leg.AddEntry(fit, ltx3Text);

    TLatex ltx4 = TLatex();
    ltx4.SetTextSize(0.035);
    char ltx4Text[256];
    snprintf(ltx4Text, sizeof(ltx4Text), "#sigma_{Core} = %.1f #pm %.1f ps", sigmaCore*1000, sigmaCoreError*1000);
    leg.AddEntry(fitCore, ltx4Text);

    TLatex ltx5 = TLatex();
    ltx5.SetTextSize(0.035);
    char ltx5Text[256];
    snprintf(ltx5Text, sizeof(ltx5Text), "#sigma_{Tail} = %.1f #pm %.1f ps", sigmaTail*1000, sigmaTailError*1000);
    leg.AddEntry(fitTail, ltx5Text);

    leg.DrawClone("Same");
	

	/* write to file */

    ofstream outfile(outputFileTxt);
    
	outfile << "fit.Entries = " << entries << ";" << endl;
	outfile << "fit.StdDev = " << stddev << ";" << endl;
	outfile << "fit.ChiSquare = " << chi2 << ";" << endl;
	outfile << "fit.NDF = " << ndf << ";" << endl;
	outfile << "fit.mean = " << mean << ";" << endl;
	outfile << "fit.mean_error = " << meanError << ";" << endl;
	outfile << "fit.sigma_core = " << sigmaCore << ";" << endl;
	outfile << "fit.sigma_core_error = " << sigmaCoreError << ";" << endl;
	outfile << "fit.sigma_tail = " << sigmaTail << ";" << endl;
	outfile << "fit.sigma_tail_error = " << sigmaTailError << ";" << endl;
	outfile << "fit.mixing = " << a << ";" << endl;
	outfile << "fit.mixing_error = " << aError << ";" << endl;
	outfile << "fit.scale = " << scale << ";" << endl;
	outfile << "fit.scale_error = " << scaleError << ";" << endl;
	outfile << "fit.sigma_all = " << sigma << ";" << endl;
	outfile << "fit.sigma_all_error = " << sigmaError << ";" << endl;
	outfile << "fit.g1_min = " << g1_min << ";" << endl;
	outfile << "fit.g1_max = " << g1_max << ";" << endl;
	outfile << "fit.g1_scale = " << par[0] << ";" << endl;
	outfile << "fit.g1_mean = " << par[1] << ";" << endl;
	outfile << "fit.g1_sigma = " << par[2] << ";" << endl;
	outfile << "fit.g1_scale_error = " << g1->GetParError(0) << ";" << endl;
	outfile << "fit.g1_mean_error = " << g1->GetParError(1) << ";" << endl;
	outfile << "fit.g1_sigma_error = " << g1->GetParError(2) << ";" << endl;

    outfile.close();
   
        
    // Print canvas to PDF with the constructed filename
    c1->Print(outputFilePdf.c_str());
    
    delete fitCore;
    delete fit;
    delete c1;
    delete h1;
    delete g1;
    

}

int double_gauss_runlist() {
    
    string inputDirectory = "Input"; // Input directory
    string outputDirectory = "ResultsROOT"; // Output directory
    
    // List of filenames to analyze
    vector<string> filenames = {
    //Photocathodes studies
"Run160-MM_metal_10 mm-120 um-PC_B4C_5_nm_on_Ti_2p4_nm-Sealed_990_mbar-A275V_C440V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run161-MM_metal_10 mm-120 um-PC_B4C_5_nm_on_Ti_2p4_nm-Sealed_990_mbar-A275V_C460V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run170-MM_metal_10 mm-120 um-PC_B4C_5_nm_on_Ti_2p4_nm-Sealed_990_mbar-A275V_C480V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run171-MM_metal_10 mm-120 um-PC_B4C_5_nm_on_Ti_2p4_nm-Sealed_990_mbar-A275V_C470V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run172-MM_metal_10 mm-120 um-PC_B4C_5_nm_on_Ti_2p4_nm-Sealed_990_mbar-A275V_C450V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run173-MM_metal_10 mm-120 um-PC_B4C_5_nm_on_Ti_2p4_nm-Sealed_990_mbar-A275V_C480V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run193-MM_metal_10 mm-120 um-PC_B4C_5_nm_no_Ti-Sealed_990_mbar-A275V_C490V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run194-MM_metal_10 mm-120 um-PC_B4C_5_nm_no_Ti-Sealed_990_mbar-A275V_C480V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run202-MM_metal_10 mm-120 um-PC_B4C_5_nm_no_Ti-Sealed_990_mbar-A275V_C470V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run205-MM_metal_10 mm-120 um-PC_B4C_5_nm_no_Ti-Sealed_990_mbar-A275V_C500V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run215-MM_metal_10 mm-120 um-PC_B4C_9_nm_on_Ti_2p4_nm-Sealed_990_mbar-A275V_C500V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run216-MM_metal_10 mm-120 um-PC_B4C_9_nm_on_Ti_2p4_nm-Sealed_990_mbar-A275V_C530V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run217-MM_metal_10 mm-120 um-PC_B4C_9_nm_on_Ti_2p4_nm-Sealed_990_mbar-A275V_C480V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run231-MM_metal_10 mm-120 um-PC_B4C_9_nm_no_Ti-Sealed_990_mbar-A275V_C480V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run232-MM_metal_10 mm-120 um-PC_B4C_9_nm_no_Ti-Sealed_990_mbar-A275V_C510V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run233-MM_metal_10 mm-120 um-PC_B4C_9_nm_no_Ti-Sealed_990_mbar-A275V_C530V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run234-MM_metal_10 mm-120 um-PC_B4C_9_nm_no_Ti-Sealed_990_mbar-A275V_C545V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run271-MM_metal_10 mm-120 um-PC_ESS_B4C_9_nm_on_Cr_2p4_nm-Sealed_990_mbar-A275V_C480V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run273-MM_metal_10 mm-120 um-PC_ESS_B4C_9_nm_on_Cr_2p4_nm-Sealed_990_mbar-A275V_C500V-MCP_cut_150_mV_and_geo_cut-SATlist.txt"
};

    int nbins = 100; // Number of bins for histogram

    // Analyze each file in the list
    for (const auto& filename : filenames) {
        cout << "Analyzing file: " << filename << endl; // Debugging print
        double_gauss(inputDirectory, outputDirectory, filename, nbins);
        cout << "Finished analyzing file: " << filename << endl; // Debugging print
    }

    gSystem->Exit(0);
    
}