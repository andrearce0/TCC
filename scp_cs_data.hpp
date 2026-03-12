#ifndef SCP_CS_DATA_HPP
#define SCP_CS_DATA_HPP

#include <iostream>
#include <vector>
#include <set>
#include <string>

// Definição da Estrutura da Instância
struct SCPCSInstance {
    int num_elementos;
    int num_subconjuntos;
    std::vector<int> custos;
    
    // Matriz de Incidência (Direta): subconjunto -> {elementos que ele cobre}
    std::vector<std::set<int>> matriz_incidencia;
    
    // Lista de Incidência (Inversa): elemento -> {subconjuntos que o cobrem}
    std::vector<std::vector<int>> lista_incidencia;
    
    // Matriz de Penalidades de Conflito (n x n)
    std::vector<std::vector<double>> matriz_conflitos; 
    int conflict_threshold; 
};
bool ler_instancia_scpcs(const std::string& nome_arquivo, SCPCSInstance& instancia, int k);
void calcular_custos_conflito(SCPCSInstance& instancia, int k);
#endif // SCP_CS_DATA_H