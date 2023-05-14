#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

void JCObinning2(){

ofstream recoedge;
ofstream recobins;
ofstream genedge;
ofstream genbins;

recoedge.open("recoedge.txt");
recobins.open("recobins.txt");
genedge.open("genedge.txt");
genbins.open("genbins.txt");

char histname[100];
double Threshold = 0.5;

int ndef = 3;
int njet = 2;
int nkappa = 10;
int njetptmn = 10;

TFile *inputData = new TFile("/home/soumyadip/Package/TUnfold/JetCharge/Prototype/Input/02May2023/PY8_bin.root");
//TFile *inputData = new TFile("/home/soumyadip/Workspace/Zigzag/JCObinning/10May2023/PY8_flat.root");

TH2D *RM[njetptmn][ndef][njet][nkappa];
TH2* hist;

genbins<<"{{{";
genedge<<"{{";
recobins<<"{{{";
recoedge<<"{{";

for(int ipt = 0; ipt < njetptmn; ipt++) {
for(int id = 0; id < ndef; id++) {
    for (int ij = 0; ij < njet; ij++) {
        for (int ik = 0; ik < nkappa; ik++) {
		if(ipt==6){
                sprintf(histname, "analyzeBasicPat/RM_jc_pt%d_eta0_d%d_j%d_k%d", ipt, id, ij, ik);
                RM[ipt][id][ij][ik] = (TH2D*) inputData->Get(histname);
                //cout<<"RM : "<<histname<<endl;

                hist = (TH2D*) RM[ipt][id][ij][ik]->Clone(); hist->RebinX(2);

                //cout <<"Definition : "<<id<<" Jet No. : "<< ij<< " Kappa : "<<ik<<" pt : "<<ipt <<endl;
	
		/*	
		cout<<"Original NBins : "<<hist->GetNbinsX()<<"  bin-edges : "<<"{";
                for (int i = 1; i <= hist->GetNbinsX()+1; i++) {
                	cout << hist->GetXaxis()->GetBinLowEdge(i)<<",";
                }
		cout<<"}";
                cout<<endl;
		*/

		double genEntriesTotal = 0.0;
		double recoEntriesTotal = 0.0;
		double purity =0.0;
		double stability =0.0;
		double array[hist->GetNbinsX()];
		int ipurity =0;
		vector<double> values;
		vector<double> reco;
		vector<double> gen;
		//cout << "Gen bin-edges : ";
		
		for(int binRec=1; binRec<= hist->GetNbinsX(); binRec++) {
			if(purity>Threshold && stability>Threshold){
				genEntriesTotal = 0.0;
				recoEntriesTotal = 0.0;
				purity = 0.0;
				stability = 0.0;
			}

			for(int binGen=1; binGen<=hist->GetNbinsX(); binGen++) {
				genEntriesTotal += hist->GetBinContent(binRec,binGen);
			}	

			for(int ireco=1; ireco<=hist->GetNbinsX(); ireco++) {
                                recoEntriesTotal += hist->GetBinContent(ireco,binRec);
                        }

			if(genEntriesTotal>0.0){
				purity += hist->GetBinContent(binRec,binRec)/genEntriesTotal;
			}

			if(recoEntriesTotal>0.0){
                                stability += hist->GetBinContent(binRec,binRec)/recoEntriesTotal;
                        }

				if(purity>Threshold && stability>Threshold){
					ipurity++;
					//cout<<"purity : "<<purity<<" stability : "<<stability<<endl;
				}
				array[binRec]=ipurity;
				
				if(binRec==1){
					genedge<<"{"<<fixed<<setprecision(2)<<hist->GetXaxis()->GetBinLowEdge(binRec)<<","; 
					values.push_back(hist->GetXaxis()->GetBinLowEdge(binRec));
				}
				if (binRec!=1 && (purity>Threshold && stability>Threshold)){
					genedge<<fixed<<setprecision(2)<<hist->GetXaxis()->GetBinLowEdge(binRec)<<","; 
					values.push_back(hist->GetXaxis()->GetBinLowEdge(binRec));
				}
				if (binRec==hist->GetNbinsX()){
					genedge<<fixed<<setprecision(2)<<hist->GetXaxis()->GetBinLowEdge(binRec+1)<<"}"; 
					values.push_back(hist->GetXaxis()->GetBinLowEdge(binRec+1));
				}
			}
			genedge<<","<<endl;
			//cout<<"Gen NBins : "<<array[hist->GetNbinsX()]+2<<endl;
			//cout<<"Reco NBins : "<<2*(array[hist->GetNbinsX()]+1)<<endl;
			//genbins<<array[hist->GetNbinsX()]+2<<",";
			//recobins<<2*(array[hist->GetNbinsX()]+2)<<",";
			//cout<<array[hist->GetNbinsX()]<<endl;

			recoedge<<"{";
			//cout<<"{";

			//Gen
			cout<<"Gen bin-edge : "<<endl;
			cout<<"{";
			for(int i=0; i<values.size()-1; i++){
				recoedge<<fixed<<setprecision(2)<<values[i]<<","<<(values[i]+values[i+1])/2<<",";
				cout<<values[i]<<",";
			}
			if(!values.empty()){cout<<values.back();}
			cout<<"}"<<endl;
			//cout<<"nGenbins "<<values.size()+1<<endl;
			genbins<<values.size()-1<<",";

			if(!values.empty()){recoedge<<values.back();}
			recoedge<<"},"<<endl;

			//Reco
			//vector<double> reco;
			cout<<"Reco bin-edge : "<<endl;
                        cout<<"{";
			for(int i=0; i<values.size()-1; i++){
				cout<<fixed<<setprecision(2)<<values[i]<<","<<(values[i]+values[i+1])/2<<",";
				reco.push_back(values[i]);
				reco.push_back((values[i]+values[i+1])/2);
			}
			if(!values.empty()){reco.push_back(values.back());cout<<values.back();}
			cout<<"}"<<endl;
			//cout<<"nRecobins "<<reco.size()+1<<endl;
			recobins<<reco.size()-1<<",";
					}
	    			}
			}
    		}
	}
genbins<<"}}}"<<endl;
genedge<<"}}"<<endl;
recobins<<"}}}"<<endl;
recoedge<<"}}"<<endl;
}
