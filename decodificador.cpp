#include "decodificador.hpp"
#include "scp_cs_data.hpp"
#include <set>
#include <vector>
#include <algorithm>
#include <limits>
#include <numeric>
#include <iterator>
#include <cmath>
#include <iostream>

using namespace std;

// Calcula pesos baseados na raridade
std::vector<double> calcular_pesos_locais(const SCPCSInstance& instancia) {
    std::vector<int> frequencia(instancia.num_elementos, 0);

    // Verifica se matriz_incidencia está preenchida
    bool matriz_ok = !instancia.matriz_incidencia.empty() && !instancia.matriz_incidencia[0].empty();

    // Se a matriz estiver vazia, tentamos usar a lista_incidencia
    if (!matriz_ok && !instancia.lista_incidencia.empty()) {
        for(int i=0; i<instancia.num_elementos; ++i) {
            frequencia[i] = instancia.lista_incidencia[i].size();
        }
    } else {
        // Fluxo normal via matriz_incidencia
        for(int j = 0; j < instancia.num_subconjuntos; j++) {
            for(int elemento : instancia.matriz_incidencia[j]) {
                if (elemento >= 0 && elemento < instancia.num_elementos) // Protecao de indice
                    frequencia[elemento]++;
            }
        }
    }

    std::vector<double> pesos(instancia.num_elementos);
    for(int i = 0; i < instancia.num_elementos; i++) {
        if (frequencia[i] > 0)
            pesos[i] = 1.0 / (double)frequencia[i];
        else
            pesos[i] = 10000.0; // Elemento nao coberto por nenhum subconjunto (erro)
    }
    return pesos;
}

double calcular_custo_solucao(const std::set<int>& subconjuntos_selecionados, const SCPCSInstance& instancia){
    double custo_total = 0.0;
    std::vector<int> selecionados(subconjuntos_selecionados.begin(), subconjuntos_selecionados.end());

    for (int j : selecionados) {
        custo_total += instancia.custos[j];
    }

    for (size_t i = 0; i < selecionados.size(); ++i) {
        int sub_i = selecionados[i];
        for (size_t j = i + 1; j < selecionados.size(); ++j) {
            int sub_j = selecionados[j];
            custo_total += instancia.matriz_conflitos[sub_i][sub_j];
        }
    }
    return custo_total;
}

double decodificar(std::vector<float> genes, 
    const SCPCSInstance& instancia, 
    std::set<int>* solucao_saida,
    const std::vector<double>& pesos_elementos) {
    int m = instancia.num_elementos;
    int n = instancia.num_subconjuntos;

    // Verifica os pesos dos elementos

    std::vector<bool> elementos_cobertos_mask(m, false);
    int elementos_cobertos_count = 0;
    std::set<int> subconjuntos_selecionados;
    std::vector<bool> ja_processado(n, false);
    double custo_total_acumulado = 0.0;
    std::vector<double> conflito_acumulado_pendente(n, 0.0);

    // Ordenar por genes
    std::vector<std::pair<float, int>> gene_prioridades;
    gene_prioridades.reserve(n);
    for (int j = 0; j < n; ++j) {
        gene_prioridades.push_back({genes[j], j});
    }
    std::sort(gene_prioridades.begin(), gene_prioridades.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    // Fase 1: Construcao Gulosa da solucao
    while (elementos_cobertos_count < m) {
        double melhor_metrica = std::numeric_limits<double>::max();
        int melhor_indice = -1;
        double custo_efetivo_do_melhor = 0.0;

        // LCR Dinamica: Se faltam poucos elementos, verifica todos os candidatos.
        // Se faltam muitos, verifique somente os X% de maior prioridade
        int range_busca = (elementos_cobertos_count > m * 0.9) ? n : std::max(1, (int)(n * 0.4));

        for (int i = 0; i < n; ++i) {
            int j = gene_prioridades[i].second;
            if (ja_processado[j]) continue;

            // Se achou um candidato e chegou ao final da lista, interrompe a busca
            if (i >= range_busca && melhor_indice != -1) break;

            // Verifica beneficio (quantos nao cobertos ele cobre + peso deles)
            double beneficio_ponderado = 0.0;
            int contagem_novos = 0;

            for (int e : instancia.matriz_incidencia[j]) {
                if (!elementos_cobertos_mask[e]) {
                    beneficio_ponderado += pesos_elementos[e];
                    contagem_novos++;
                }
            }

            if (contagem_novos == 0) continue;

            double penalidade_conf = conflito_acumulado_pendente[j];
            double custo_efetivo = (double)instancia.custos[j] + penalidade_conf;

            // Heuristica
            if (beneficio_ponderado < 1e-9) beneficio_ponderado = 1e-9;
            double metrica_gulosa = custo_efetivo / beneficio_ponderado;

            if (metrica_gulosa < melhor_metrica) {
                melhor_metrica = metrica_gulosa;
                melhor_indice = j;
                custo_efetivo_do_melhor = custo_efetivo;
            }
        }

        // Se nao achou ninguem nesta rodada, encerra a Fase 1
        if (melhor_indice == -1) break;

        // Aplica escolha
        custo_total_acumulado += custo_efetivo_do_melhor;
        subconjuntos_selecionados.insert(melhor_indice);
        ja_processado[melhor_indice] = true;

        for (int e : instancia.matriz_incidencia[melhor_indice]) {
            if (!elementos_cobertos_mask[e]) {
                elementos_cobertos_mask[e] = true;
                elementos_cobertos_count++;
            }
        }

        // Atualiza penalidades
        for (int k = 0; k < n; ++k) {
            if (!ja_processado[k]) {
                conflito_acumulado_pendente[k] += instancia.matriz_conflitos[melhor_indice][k];
            }
        }
    }

    //FASE 2: REPARO
    //Se a fase gulosa terminou mas existem elementos descobertos
   if (elementos_cobertos_count < m) {
        // Itera sobre cada elemento para ver se foi coberto
        for (int i = 0; i < m; ++i) {
            if (!elementos_cobertos_mask[i]) {
                int melhor_candidato_emergencia = -1;
                double menor_custo_emergencia = std::numeric_limits<double>::max();

                // Verifica quem cobre o elemento i
                const auto& cobridores = instancia.lista_incidencia[i];

                if (cobridores.empty()) {
                    // Se ninguem cobre 'i', a instância é inválida. Retorna Infinito.
                    return std::numeric_limits<double>::infinity();
                }

                // Escolhe o mais barato dentre os que cobrem 'i'
                for (int cand : cobridores) {
                    if (ja_processado[cand]) continue; // Já está na solução (bug lógico se entrar aqui, mas segurança)

                    double custo_cand = instancia.custos[cand] + conflito_acumulado_pendente[cand];
                    if (custo_cand < menor_custo_emergencia) {
                        menor_custo_emergencia = custo_cand;
                        melhor_candidato_emergencia = cand;
                    }
                }

                // Se achou candidato, adiciona a solucao
                if (melhor_candidato_emergencia != -1) {
                    subconjuntos_selecionados.insert(melhor_candidato_emergencia);
                    ja_processado[melhor_candidato_emergencia] = true;
                    custo_total_acumulado += menor_custo_emergencia;

                    // Atualiza coberturas
                    for (int e : instancia.matriz_incidencia[melhor_candidato_emergencia]) {
                        if (!elementos_cobertos_mask[e]) {
                            elementos_cobertos_mask[e] = true;
                            elementos_cobertos_count++;
                        }
                    }

                    // Atualiza conflitos
                    for (int k = 0; k < n; ++k) {
                        if (!ja_processado[k]) {
                            conflito_acumulado_pendente[k] += instancia.matriz_conflitos[melhor_candidato_emergencia][k];
                        }
                    }
                }
            }
        }
    }

    if(solucao_saida != nullptr)
        *solucao_saida = subconjuntos_selecionados;

    //Retorna custo real
    return calcular_custo_solucao(subconjuntos_selecionados, instancia);
}
