#include <stdio.h>
#include <stdlib.h>
#include <TMath.h>
#include <vector>
#include <string>
#include <TTree.h>
#include <TFile.h>
#include "TCanvas.h"

using namespace std;

//Script para visualizar o traco de um evento

void TraceView (int runnumber, int IdEvento) //numero do run, numero do evento
{
    
    string nomearq = Form("/home/barbaras/red-deconv/barb/trees_v2/mytreefile_%d.root", runnumber); //Atualize o caminho para o seu pc!!
    TFile *file = new TFile(nomearq.c_str(),"READ"); //Abrindo o arquivo em modo de leitura 

    //Checando se o arquivo escolhido realmente existe
    if ( file->IsZombie() ) 
    {   printf("RedFile does not exist \n"); 
        return 0; 
    }

    //Carregando a arvore para uma variavel do script
    TTree *trace = (TTree*)file->Get("treetraces"); 

    //Criando um vetor para armazenar os valores de energia
    vector<float> *energy = nullptr; 
    
    //Associando o galho da arvore ao vetor criado
    trace->SetBranchAddress("energy",&energy); 


    if (IdEvento > trace->GetEntries()) //GetEntries() -> retorna o número total de entradas na árvore, a.k.a. num de eventos
    {
        cout<<"Erro. esse evento nao existe"<<endl; //Checando se o evento existe!
        return 0;
    }

    //Carregando as informacoes do evento para a arvore
    trace->GetEntry(IdEvento); 

    //Criando histograma e grafico para desenhar o traco 
    TH1F *histo = new TH1F("histo","",100, 0., 100.);
    TGraph *graph = new TGraph(energy->size());


    
    //Preenchendo o histograma e o grafico com o traco do evento
     for (int i=0; i<(int) energy-> size(); i++) 
    {
        histo->Fill(i,energy->at(i)); //at(i) -> função que pega a entrada do vetor na posicao i
        graph->SetPoint(i+1, i, energy->at(i));
    }

    trace->ResetBranchAddresses(); //Resetando memorias da arvore


    //Desenhando histograma e grafico 

    auto c1 = new TCanvas("c1","c1",200,10,700,500);
    histo->Draw();


    auto c2 = new TCanvas("c2","c2",200,10,700,500);
    graph->Draw("ALP"); //a = axis, L = linha ligando os pontos, p = pontos

}