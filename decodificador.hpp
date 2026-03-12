#ifndef DECODIFICADOR_HPP
#define DECODIFICADOR_HPP

#include "scp_cs_data.hpp"
#include <vector>

double decodificar(std::vector<float> genes, const SCPCSInstance& instancia, std::set<int>* solucao_saida = nullptr, const std::vector<double>& pesos_elementos = {});
std::vector<double> calcular_pesos_locais(const SCPCSInstance& instancia);
#endif // DECODIFICADOR_HPP
