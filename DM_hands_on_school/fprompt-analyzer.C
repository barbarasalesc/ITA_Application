#include <TMath.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <string>
#include <TTree.h>

// Função principal que calcula o Fprompt recebendo as janelas de integração curta (wshort) e longa (wlong)
void FpromptCalc(int wshort, int wlong) {
    
    // Cria um histograma 2D para plotar o Fprompt em função da Energia
    // Eixos: 1000 bins de energia (0 a 10000) x 110 bins de Fprompt (0.4 a 1.0)
    TH2D *hFpromptE = new TH2D("hFpromptE","2D Histogram",1000,0.,10000.,110,0.4,1.);
    
    // Abre um arquivo de texto contendo a lista de "runs" (corridas de tomada de dados) a serem analisados
    std::ifstream f("analysis/selection.runs",std::ifstream::in);

    // Loop principal sobre as corridas listadas no arquivo de texto
    while (true)
    {
        int run;
        f >> run; // Lê o número do run
        if ( f.eof() ) break; // Sai do loop quando chegar ao fim do arquivo

        // Abre o arquivo ROOT correspondente à corrida atual
        TFile *file = new TFile(Form("trees_v2/mytreefile_%d.root",run), "READ");
        
        // Verifica se o arquivo está corrompido ou não existe para evitar falhas de execução
        if ( file->IsZombie() ) { printf("RedFile does not exist \n"); continue; }
        
        // Carrega a árvore de dados (TTree) que contém os traços (traces) do evento
        TTree *trace = (TTree*)file->Get("treetraces");
        
        // Inicialização de ponteiros e variáveis para ler os branches (ramos) da TTree
        std::vector<float> *energy  = nullptr;
        int sa_start, sa_end;
        
        // Associa as variáveis locais aos branches da árvore
        trace->SetBranchAddress("energy",&energy);
        trace->SetBranchAddress("sa_start",&sa_start);
        trace->SetBranchAddress("sa_end",&sa_end);

        // Loop sobre todos os eventos (entradas) dentro da árvore da corrida atual
        for (int ientry=0; ientry<trace->GetEntries(); ientry++)
        {
            trace->GetEntry(ientry); // Carrega os dados do evento atual para a memória
            
            // Imprime no terminal o progresso a cada 1000 eventos processados
            if ( ientry%1000==0 ) printf("Run %d - Entry %d \n",run,ientry);

            // Variáveis para acumular a carga nas janelas de tempo e calcular o Fprompt
            double qshort = 0;
            double qlong = 0;
            double totalenergy = 0;
            double energytoplot = 0;
            double fprompt = 0;

            // Integra a carga no componente rápido (prompt) do sinal
            for (int sa=0; sa<wshort; sa++)
            {
                qshort += energy->at(sa);
            } 
            // Integra a carga no componente longo do sinal
            for (int sa=0; sa<wlong; sa++)
            {
                qlong += energy->at(sa); 
            } 
            // Integra a carga total de todo o traço gravado
            for (int sa=0; sa<(int) energy->size(); sa++)
            {
                totalenergy += energy->at(sa);
            } 
            
            // Calcula o parâmetro Fprompt (razão entre a carga rápida e a carga total/longa)
            fprompt = qshort/qlong;
            
            // Inverte o sinal da energia total para plotagem (presumindo sinais de polaridade negativa no detector)
            energytoplot = -totalenergy;
            
            // Preenche o histograma 2D com a energia do evento e o seu Fprompt calculado
            hFpromptE->Fill(energytoplot, fprompt);
        }
        
        // Limpa os endereços dos branches antes de carregar a árvore da próxima corrida
        trace->ResetBranchAddresses();
    }

    // Cria um Canvas de alta resolução e desenha o histograma 2D com mapa de cores (COLZ)
    TCanvas *c1 = new TCanvas ("c1","c1",1920,1080);
    hFpromptE->Draw("COLZ");

    // Faz a projeção do histograma 2D no eixo Y (obtendo a distribuição 1D apenas do parâmetro Fprompt)
    TH1D *hfprompt = (TH1D*) hFpromptE->ProjectionY("hfprompt");
    
    // Cria um segundo Canvas e desenha a projeção 1D
    TCanvas *c2 = new TCanvas ("c2","c2",1920,1080);
    hfprompt->Draw();
}
