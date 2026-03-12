#include "scp_cs_data.hpp"
#include "decodificador.hpp"
#include <iostream>
#include <string>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <future>
#include <chrono>
#include <set>
#include <limits>
#include <cmath>
#include <thread>

using namespace std;

struct BRKGAConfig {
    int tamanho_populacao;
    int tamanho_elite;
    int tamanho_mutantes;
    int num_ilhas;
    int freq_migracao;
    int num_migrantes;
    double rho;
    double tempo_limite;
};

struct Cromossomo {
    std::vector<float> genes;
    double fitness = 0;
};

// Geradores de numeros aleatorios
static std::mt19937_64& get_rng() {
    static std::mt19937_64 rng(std::random_device{}());
    return rng;
}

float gerar_numero_aleatorio() {
    std::uniform_real_distribution<double> dist(0.01, 0.99);
    return dist(get_rng());
}

int gerar_indice_aleatorio(int max_exclusive) {
    if (max_exclusive <= 0) return 0;
    std::uniform_int_distribution<int> dist_int(0, max_exclusive - 1);
    return dist_int(get_rng());
}

// Parametros para o brkga
BRKGAConfig configurar_brkga_automaticamente(const SCPCSInstance& instancia, double tempo_limite_arg) {
    BRKGAConfig cfg;
    int n = instancia.num_subconjuntos;

    // Numero de ilhas
    unsigned int nucleos = std::thread::hardware_concurrency();
    if (nucleos == 0) nucleos = 4;

    if (n < 100) {
        cfg.num_ilhas = 2;
    } else {
        cfg.num_ilhas = std::min((int)nucleos, 4);
        if (cfg.num_ilhas < 2) cfg.num_ilhas = 2;
    }

    // Tamanho da Populacao
    int pop_total = (int)n * 1.5;
    if (pop_total > 1800) pop_total = 1800; // Pop maxima = 1800 cromossomos

    while (pop_total % cfg.num_ilhas != 0) pop_total++;

    cfg.tamanho_populacao = pop_total;

    // Tamanho da Elite e Mutantes
    cfg.tamanho_elite = (int)(0.15 * cfg.tamanho_populacao);
    cfg.tamanho_mutantes = (int)(0.25 * cfg.tamanho_populacao); 

    while (cfg.tamanho_elite % cfg.num_ilhas != 0) cfg.tamanho_elite++;

    cfg.rho = 0.60;
    cfg.freq_migracao = 20;
    cfg.num_migrantes = (cfg.num_ilhas > 1) ? 10 : 0;
    cfg.tempo_limite = tempo_limite_arg;

    return cfg;
}

double calcular_fitness_real(const std::set<int>& solucao, const SCPCSInstance& instancia) {
    double custo_total = 0.0;
    for (int s : solucao) custo_total += instancia.custos[s];

    if (!solucao.empty()) {
        std::vector<int> vec_sol(solucao.begin(), solucao.end());
        for (size_t i = 0; i < vec_sol.size(); ++i) {
            for (size_t j = i + 1; j < vec_sol.size(); ++j) {
                int u = vec_sol[i];
                int v = vec_sol[j];
                custo_total += instancia.matriz_conflitos[u][v];
            }
        }
    }
    return custo_total;
}

void aplicar_poda(std::set<int>& solucao, const SCPCSInstance& instancia) {
    if (solucao.empty()) return;
    std::vector<int> contagem_cobertura(instancia.num_elementos, 0);

    for (int subconjunto : solucao) {
        for (int elemento : instancia.matriz_incidencia[subconjunto]) {
            contagem_cobertura[elemento]++;
        }
    }

    std::vector<int> candidatos(solucao.begin(), solucao.end());
    std::sort(candidatos.begin(), candidatos.end(), [&](int a, int b) {
        return instancia.custos[a] > instancia.custos[b];
    });

    for (int subconjunto_teste : candidatos) {
        bool essencial = false;
        for (int elemento : instancia.matriz_incidencia[subconjunto_teste]) {
            if (contagem_cobertura[elemento] <= 1) {
                essencial = true;
                break;
            }
        }

        if (!essencial) {
            solucao.erase(subconjunto_teste);
            for (int elemento : instancia.matriz_incidencia[subconjunto_teste]) {
                contagem_cobertura[elemento]--;
            }
        }
    }
}

double decodificar_com_poda(const std::vector<float>& genes, const SCPCSInstance& instancia, const std::vector<double>& pesos) {
    std::set<int> solucao;
    double fit_guloso = decodificar(genes, instancia, &solucao, pesos);

    if (fit_guloso == std::numeric_limits<double>::infinity()) {
        return fit_guloso;
    }

    aplicar_poda(solucao, instancia);
    return calcular_fitness_real(solucao, instancia);
}

// GESTAO DE POPULAÇÃO
vector<Cromossomo> gerar_populacao_inicial(const SCPCSInstance& instancia, int tamanho_populacao) {
    vector<Cromossomo> populacao;
    populacao.reserve(tamanho_populacao);

    for (int i = 0; i < tamanho_populacao; i++) {
        Cromossomo novo;
        novo.genes.resize(instancia.num_subconjuntos);
        for (int j = 0; j < instancia.num_subconjuntos; j++) {
            novo.genes[j] = gerar_numero_aleatorio();
        }
        populacao.push_back(novo);
    }
    return populacao;
}

void aplicar_fitness_paralela(vector<Cromossomo>& populacao, int indice_inicio_novos, const SCPCSInstance& instancia, const std::vector<double>& pesos) {
    std::vector<std::future<double>> futuros;
    int num_novos = populacao.size() - indice_inicio_novos;
    if(num_novos <= 0) return;
    futuros.reserve(num_novos);

    for (size_t i = indice_inicio_novos; i < populacao.size(); ++i) {
        futuros.push_back(std::async(std::launch::async, decodificar_com_poda, populacao[i].genes, std::cref(instancia), std::cref(pesos)));
    }

    int fut_idx = 0;
    for (size_t i = indice_inicio_novos; i < populacao.size(); ++i) {
        populacao[i].fitness = futuros[fut_idx].get();
        fut_idx++;
    }
}

double brkga(SCPCSInstance& instancia, BRKGAConfig config) {

    // Configuração base das ilhas
    int tam_pop_ilha = config.tamanho_populacao / config.num_ilhas;
    int tam_elite_ilha = config.tamanho_elite / config.num_ilhas;
    int tam_mutantes_ilha = config.tamanho_mutantes / config.num_ilhas;

    cout << "tamanho elite ilha: " << tam_elite_ilha << endl;
    cout << "tamanho mutantes ilha: " << tam_mutantes_ilha << endl;

    // Proteções
    if (tam_elite_ilha < 1) tam_elite_ilha = 1;
    if (tam_pop_ilha < tam_elite_ilha + 5) tam_pop_ilha = tam_elite_ilha + 10;

    std::vector<double> pesos = calcular_pesos_locais(instancia);

    // Inicialização das ilhas
    vector<vector<Cromossomo>> ilhas(config.num_ilhas);
    for(int k=0; k < config.num_ilhas; k++) {
        ilhas[k] = gerar_populacao_inicial(instancia, tam_pop_ilha);
        aplicar_fitness_paralela(ilhas[k], 0, instancia, pesos);
        std::sort(ilhas[k].begin(), ilhas[k].end(), [](const Cromossomo& a, const Cromossomo& b) {
            return a.fitness < b.fitness;
        });
    }

    int cont_geracao = 0;
    double melhor_solucao_global = std::numeric_limits<double>::max();
    Cromossomo melhor_cromossomo_global;

    // Melhor inicial
    for(int k=0; k < config.num_ilhas; k++){
        if(ilhas[k][0].fitness < melhor_solucao_global){
            melhor_solucao_global = ilhas[k][0].fitness;
            melhor_cromossomo_global = ilhas[k][0];
        }
    }

    int geracao_melhor_solucao = 0;
    double tempo_melhor_solucao = 0.0;
    auto inicio = std::chrono::steady_clock::now();

    while (true) {
        auto tempo_atual = std::chrono::steady_clock::now();
        double tempo_total = std::chrono::duration_cast<std::chrono::duration<double>>(tempo_atual - inicio).count();

        // Loop por ilha
        for(int k=0; k < config.num_ilhas; k++) {

            vector<Cromossomo>& pop_atual = ilhas[k];
            vector<Cromossomo> nova_populacao;
            nova_populacao.reserve(tam_pop_ilha);

            // Copia Elite para prox populacao
            for(int i=0; i<tam_elite_ilha; ++i)
                nova_populacao.push_back(pop_atual[i]);

            // Gera Mutantes
            vector<Cromossomo> mutantes = gerar_populacao_inicial(instancia, tam_mutantes_ilha);
            nova_populacao.insert(nova_populacao.end(), mutantes.begin(), mutantes.end());

            // Crossover
            // Escolhe pai1 elite
            int vagas_restantes = tam_pop_ilha - nova_populacao.size();
            for (int i = 0; i < vagas_restantes; i++) {
                const Cromossomo& pai1 = pop_atual[gerar_indice_aleatorio(tam_elite_ilha)];

                // Escolhe pai2 nao-elite
                int idx_pai2 = tam_elite_ilha + gerar_indice_aleatorio(tam_pop_ilha - tam_elite_ilha);
                const Cromossomo& pai2 = pop_atual[idx_pai2];

                Cromossomo filho;
                filho.genes.resize(instancia.num_subconjuntos);

                for (int j = 0; j < instancia.num_subconjuntos; j++) {
                    if (gerar_numero_aleatorio() <= config.rho) filho.genes[j] = pai1.genes[j];
                    else filho.genes[j] = pai2.genes[j];
                }
                nova_populacao.push_back(filho);
            }

            ilhas[k] = nova_populacao;
            aplicar_fitness_paralela(ilhas[k], tam_elite_ilha, instancia, pesos);

            std::sort(ilhas[k].begin(), ilhas[k].end(), [](const Cromossomo& a, const Cromossomo& b) {
                return a.fitness < b.fitness;
            });
            // Atualiza Global
            if (ilhas[k][0].fitness < melhor_solucao_global) {
                melhor_solucao_global = ilhas[k][0].fitness;
                melhor_cromossomo_global = ilhas[k][0];
                geracao_melhor_solucao = cont_geracao;
                tempo_melhor_solucao = tempo_total;

                cout << ">>> Melhoria Global (Ilha " << k << "): " << melhor_solucao_global
                     << " [Gen " << cont_geracao << "] apos " << tempo_melhor_solucao << "s" << endl;
            }
        }

        cont_geracao++;

        // Migracao entre ilhas
        if (config.num_ilhas > 1 && cont_geracao % config.freq_migracao == 0) {
            vector<Cromossomo> buffer;
            for(int m=0; m<config.num_migrantes; m++) buffer.push_back(ilhas[config.num_ilhas-1][m]);

            for(int k=0; k<config.num_ilhas; k++) {
                vector<Cromossomo> prox;
                if (k < config.num_ilhas - 1) {
                     for(int m=0; m<config.num_migrantes; m++) prox.push_back(ilhas[k][m]);
                }
                for(int m=0; m<config.num_migrantes; m++) {
                    ilhas[k][tam_pop_ilha - 1 - m] = buffer[m];
                }
                buffer = prox;
                std::sort(ilhas[k].begin(), ilhas[k].end(), [](const Cromossomo& a, const Cromossomo& b) {
                    return a.fitness < b.fitness;
                });
            }
        }

        if (tempo_total >= config.tempo_limite || cont_geracao >= 150) {
            cout << "Tempo limite ou maximo de geracoes atingido." << endl;
            cout << "Tempo total: " << tempo_total << "s" << endl;
            cout << "Geracoes: " << cont_geracao << endl;
            break;
        }

    }

    cout << "\n=== RESULTADO FINAL ===" << endl;
    cout << "Melhor Custo: " << melhor_solucao_global << endl;
    cout << "Encontrado em: " << tempo_melhor_solucao << "s (Geracao " << geracao_melhor_solucao << ")" << endl;
    cout << "Total geracoes: " << cont_geracao << endl;
    std::set<int> conjunto_solucao;
    decodificar(melhor_cromossomo_global.genes, instancia, &conjunto_solucao, pesos);
    aplicar_poda(conjunto_solucao, instancia);

    cout << "Elementos Solucao (" << conjunto_solucao.size() << "): ";
    for(auto elemento : conjunto_solucao) cout << elemento + 1 << " ";
    cout << endl << endl;

    return melhor_solucao_global;
}

int main(int argc, char* argv[]) {
    string nome_arquivo;
    double tempo_limite = 600.0;

    if (argc < 2) {
        nome_arquivo = "instancias//scp41-3.txt";
    } else {
        nome_arquivo = argv[1];
        if (argc >= 3) tempo_limite = std::stod(argv[2]);
    }

    SCPCSInstance inst;
    int k_threshold = 1;

    cout << "Carregando instancia: " << nome_arquivo << "..." << endl;
    if(!ler_instancia_scpcs(nome_arquivo, inst, k_threshold)) {
        cerr << "Erro fatal: Nao foi possivel ler o arquivo." << endl;
        return 1;
    }
    calcular_custos_conflito(inst, k_threshold);

    BRKGAConfig config = configurar_brkga_automaticamente(inst, tempo_limite);

    cout << "=== BRKGA MULTI-ILHAS (STANDARD) ===" << endl;
    cout << "Subconjuntos: " << inst.num_subconjuntos << endl;
    cout << "Ilhas:        " << config.num_ilhas << endl;
    cout << "Pop/Ilha:     " << config.tamanho_populacao / config.num_ilhas << endl;
    cout << "Tempo Limite: " << config.tempo_limite << "s" << endl;
    cout << "====================================" << endl;

    int numero_testes = 10; //Numero de vezes que a instancia sera executada
    for(int teste=1; teste <= numero_testes; teste++) {
        cout << "\n--- Teste " << teste << " de " << numero_testes << " ---" << endl;
        brkga(inst, config);
    }

    system("pause");
    return 0;
}
