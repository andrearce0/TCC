# SCP-CS-BRKGA
Implementação desenvolvida como parte do TCC de minha graduação em Ciência da Computação. Consiste em um algoritmo BRKGA implementado em C++ para resolver instâncias do Set Covering Problem with Conflicts on Subsets. O projeto consiste em 3 principais arquivos:
- 1 - scp_cs_data.cpp: Arquivo que implementa as funções de leitura do arquivo .txt que contém as informações de uma instância SCP (os arquivos utilizados encontram-se na pasta "instancias"). Uma vez lido o arquivo, as informações da instância são armazenadas e utilizadas por outros arquivos do projeto.
  
- 2 - decodificador.cpp: Implementa a função de decodificar cromossomo (torná-lo em solução viável) e a função de busca local, que é aplicada ao final do brkga no melhor indivíduo encontrando, visando remover subconjuntos redundantes.
  
- 3 - brga.cpp: O principal arquivo do projeto, que implementa o brkga de fato. Em sua função main são definidos os parâmetros do brkga, como o nome da instância a ser lida, tamanho da população, tamanho do conjunto elite, etc.

Após definir os parâmetros na função main do arquivo brkga.cpp, para gerar o .exe do projeto basta executar no terminal o seguinte comando:

- g++ brkga.cpp scp_cs_data.cpp decodificador.cpp -Iinclude -o brkga
