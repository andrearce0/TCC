#include "scp_cs_data.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <cmath>

using namespace std;

bool ler_instancia_scpcs(const std::string& nome_arquivo, SCPCSInstance& instancia, int k) {
    ifstream arquivo(nome_arquivo);
    if (!arquivo.is_open()) {
        cerr << "Erro ao abrir o arquivo: " << nome_arquivo << endl;
        return false;
    }

    instancia.conflict_threshold = k;
    string linha;
    
    // Leitura da primeira linha: m (numero de elementos) e n (numero de subconjuntos)
    if (!getline(arquivo, linha)) {
        cerr << "Erro: Arquivo vazio ou nao foi possivel ler a primeira linha (m e n)." << endl;
        return false;
    }
    stringstream ss_mn(linha);
    ss_mn >> instancia.num_elementos >> instancia.num_subconjuntos;

    int n = instancia.num_subconjuntos;
    int m = instancia.num_elementos;
    
    // Inicializa estruturas com base em n e m
    instancia.custos.resize(n); //Vetor de custos c_j
    instancia.matriz_incidencia.resize(n); //Matriz de Incidência (elementos na linha j são cobertos pelo subconjunto j)
    instancia.lista_incidencia.resize(m); //Lista de Incidência (elementos na linha i são os subconjuntos que cobrem o elemento i)

    // Cria uma string para consumir todos os custos e IDs.
    string buffer_dados;
    while (getline(arquivo, linha)) {
        buffer_dados += linha;
        buffer_dados += " "; //Espaco para separar números entre linhas
    }
    
    stringstream ss_buffer(buffer_dados);

    //Leitura dos custos (c_j)
    for (int j = 0; j < n; ++j) {
        if (!(ss_buffer >> instancia.custos[j])) {
            cerr << "Erro: Arquivo terminou ou formato incorreto durante a leitura dos " << n << " custos." << endl;
            return false;
        }
    }

    // Leitura das linhas de incidencia
    for (int i = 0; i < m; ++i) {
        int num_coberturas;
        
        // Le o número de subconjuntos que cobrem o elemento i
        if (!(ss_buffer >> num_coberturas)) {
            cerr << "Erro: Faltam dados de 'num_coberturas' para o elemento " << i + 1 << endl;
            return false;
        }

        int subconjunto_id_1;
        //percorre a lista de IDs de subconjunto (num_coberturas vezes)
        for (int k_id = 0; k_id < num_coberturas; ++k_id) {
            
            if (!(ss_buffer >> subconjunto_id_1)) {
                cerr << "Erro: Arquivo terminou ou formato incorreto durante a leitura dos IDs para o elemento " << i + 1 << endl;
                return false;
            }

            //arquivo é lido em base 1 e convertido para base 0
            int subconjunto_id_0 = subconjunto_id_1 - 1;
            int elemento_id_0 = i;
            
            //verificar se os ids estao no limite estabelecido
            if (subconjunto_id_0 < 0 || subconjunto_id_0 >= n) {
                 cerr << "Erro de índice: ID de subconjunto (" << subconjunto_id_1 << ") fora dos limites (1 a " << n << "). SEGFUALT em indice " << subconjunto_id_0 << "." << endl;
                 return false;
            }

            //construcao da matriz incidencia
            instancia.matriz_incidencia[subconjunto_id_0].insert(elemento_id_0); 
            instancia.lista_incidencia[elemento_id_0].push_back(subconjunto_id_0);
        }
    }

    std::cout << "Leitura do arquivo texto concluida com sucesso. Elementos (m): " << instancia.num_elementos 
         << ", Subconjuntos (n): " << instancia.num_subconjuntos << endl;
    
    return true;
}

void calcular_custos_conflito(SCPCSInstance& instancia, int k) {
    instancia.conflict_threshold = k;
    int n = instancia.num_subconjuntos;

    //redimensiona e inicializa a Matriz de Conflitos (todos a 0.0)
    instancia.matriz_conflitos.resize(n, vector<double>(n, 0.0));

    double max_ratio = 0.0;
    
    for (int j = 0; j < n; ++j) {
        // custos[j] é c_j, e matriz_incidencia[j].size() é a cardinalidade do subconjunto j
        int card_subconjunto = instancia.matriz_incidencia[j].size();
        
        if (card_subconjunto > 0) {
            double ratio = (double)instancia.custos[j] / card_subconjunto;
            if (ratio > max_ratio) {
                max_ratio = ratio;
            }
        }
    }
    
    //a maior relacao custo/cardinalidade representa o custo de conflito unitario
    int coeff = max(1, (int)round(max_ratio)); 
    
    //iteração sobre todos os pares de subconjuntos (i, j) para calcular os custos de conflito
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            
            const std::set<int>& sub_i = instancia.matriz_incidencia[i];
            const std::set<int>& sub_j = instancia.matriz_incidencia[j];

            //calcular quantos elementos os dois subconjuntos têm em comum
            int intersection_size = 0;
            
            //verifica quais elementos do menor conjunto estão no maior
            if (sub_i.size() < sub_j.size()) {
            for (int elemento : sub_i) {
                if (sub_j.count(elemento)) {
                    intersection_size++;
                }
            }
            } else {
                for (int elemento : sub_j) {
                if (sub_i.count(elemento)) {
                    intersection_size++;
                }
                }
            }
            
            //o valor de conflito é baseado em quantos elementos em comum excedem o limiar k
            int conflict_size = max(0, intersection_size - k);
        
            if (conflict_size > 0) {
                //o custo unitario de conflito é multiplicado pelo número de conflitos
                double conflict_cost = (double)coeff * conflict_size;
                
                //armazenar o custo na matriz
                instancia.matriz_conflitos[i][j] = conflict_cost;
                instancia.matriz_conflitos[j][i] = conflict_cost;
            }
        }
    }
}
