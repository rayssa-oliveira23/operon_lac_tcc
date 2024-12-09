/*
    Simulation of lac operon model with Gillespie's Algorithm
    Rayssa Oliveira Santos
    
    argv[1] - maximum simulation time

    Verson with binary file output

    Compilation command: 
    g++ ssa_lac_operon_iptg.cpp -o ssa_lac_operon_iptg
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <string>

void SSA(double tfinal, std::vector<std::vector<double>> V, std::vector<double> X, std::vector<double> c, double I_ex);
std::vector<double> propensity(std::vector<double> X, std::vector<double> c, double I_ex);

int main(int argc, char* argv[])
{   
    double tfinal = std::stod(argv[1]);         // Maximum simulation time
    double I_ex = std::stod(argv[2]);           // Iex volume (nM) | Number of molecules of Iex = 9635 (20 nM).

    // tfinal do Stamatakis = 10^4 - 10^6
    std::vector<double> X, c;
    std::vector<std::vector<double>> V;

    //-- Lac operon
    // Initial quantity of molecules
    //
    //  M_R,   R,  R2,    O,  R2O,  I, I2R2, M_Y,  Y, YI_ex
    X = {0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};          


    // Stoichiometric matrix
    V = {{ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0},   // M_R
         { 0,  1, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0},   // R
         { 0,  0,  1, -1, -1,  1, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0},   // R2
         { 0,  0,  0,  0, -1,  1,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},   // O
         { 0,  0,  0,  0,  1, -1,  0,  0, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},   // R2O
         { 0,  0,  0,  0,  0,  0, -2,  2, -2,  2,  0,  0,  0,  0,  0,  1,  1, -1,  0,  0,  0,  0,  0,  1,  2},   // I
         { 0,  0,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1},   // I2R2
         { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0},   // M_Y
         { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1,  1,  1,  0,  0,  0,  0,  0,  0, -1,  0,  0},   // Y 
         { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  0,  0,  0,  0,  0,  0,  0, -1,  0}};  // YI_ex
    //     1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25

    // Reaction parameters
    c = {0.23, 15., 50., 1e-3, 960., 2.4, 3e-7, 12., 3e-7, 4.8e3, 0.5, 0.01, 30., 0.12, 0.1, 6e4, 0.92, 0.92, 0.462, 0.462, 0.2, 0.2, 0.2, 0.2, 0.2}; 

    SSA(tfinal, V, X, c, I_ex);

    return 0;
}

void SSA(double tfinal, std::vector<std::vector<double>> V, std::vector<double> X, std::vector<double> c, double I_ex)
{
    int j = 0, m = c.size(), n = X.size();
    double asum, aux, rand1, rand2, tau, t=0;
    std::vector<double> a, cumsum;
    bool flag = false;

    // Gerador de números aleatórios Mersenne Twister
    std::mt19937_64 mt(time(nullptr));
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Para exportar os resultados para um arquivo
    std::ofstream xFile("x.bin", std::ios::binary | std::ios::app);
    std::ofstream tFile("t.bin", std::ios::binary | std::ios::app);

    while(t < tfinal)
    {  
        //---- Passo 1: determinar qual reação vai ocorrer (reação j)
        // Função de propensidade, um item para cada reação possível
        a = propensity(X, c, I_ex);

        // Soma dos itens de a
        asum = 0.0;
        for (double item : a) asum += item;

        // Soma acumulativa dos valores de a divididor por asum
        aux = 0;
        for (size_t i = 0; i < m; i++) {
            aux += a[i];
            cumsum.push_back(aux/asum);
        }

        // Geração do primeiro número aleatório
        rand1 = dis(mt);

        // Determinar o *índice* do menor dos itens de cumsum que são maiores que esse número aleatório
        for (size_t i = 0; i < m; i++)
        {
            if (cumsum[i] > rand1){j = i; break;}
        }
        cumsum.clear();

        //---- Passo 2: determinar quando a reação j irá ocorrer (tau)
        // Geração do segundo número aleatório
        rand2 = dis(mt);

        // O número aleatório rand2 passa por transformação inversa para que tau siga uma distribuição exponencial
        tau = log(1/rand2)/asum;
    
        //---- Passo 3: atualizar vetor de estados após ocorrer a reação j, isto é, X += coluna j da matriz estequiométrica V
        xFile.write(reinterpret_cast<const char*>(X.data()), n * sizeof(double));
        for (size_t i = 0; i < n; ++i)
        {
            X[i] += V[i][j];
        }

        tFile.write(reinterpret_cast<const char*>(&t), sizeof(t));

        //---- Passo 4: atualizar o tempo
        t += tau;
        std::cout << t << "\n";
    }
    xFile.close();
    tFile.close();
}

std::vector<double> propensity(std::vector<double> X, std::vector<double> c, double I_ex)
{   
    int m = 25;                 // Number of reactions
    double V = 8.0e-16;         // E. coli volume (L)
    double N_A = 6.0221367e14;  // Avogadro's number (nmol^-1)

    std::vector<double> a(m);

    // 1. LacI transcription                            0 -> M_R                (k_sMR)
    a[0] = V*N_A*c[0];
    
    // 2. LacI monomer translation                      M_R -> M_R + R          (k_sR)
    a[1] = c[1]*X[0];

    // 3. LacI dimerization                             2R -> R2                (k_2R)
    a[2] = ((c[2])/(V*N_A))*X[1]*(X[1]-1);

    // 4. LacI dimer dissociation                       R2 -> 2R                (k_-2R)
    a[3] = c[3]*X[2];

    // 5. Association (repression)                      R2 + O -> R2O           (k_r)
    a[4] = ((c[4])/(V*N_A))*X[2]*X[3];
    
    // 6. Dissociation (repression)                     R2O -> R2 + O           (k_-r)
    a[5] = c[5]*X[4];

    // 7. 1st derepression mechanism (association)      2I + R2 -> I2R2         (k_dr1)
    a[6] = ((c[6])/(V*N_A))*X[2]*X[5]*(X[5]-1);

    // 8. 1st derepression mechanism (dissociation)     I2R2 -> 2I + R2         (k_-dr1)
    a[7] = c[7]*X[6];

    // 9. 2st derepression mechanism (association)      2I + R2O -> I2R2 + O    (k_dr2)
    a[8] = ((c[8])/(V*N_A))*X[4]*X[5]*(X[5]-1);
 
    // 10. 2st derepression mechanism (dissociation)    I2R2 + O -> 2I + R2O    (k_-dr2)
    a[9] = ((c[9])/(V*N_A))*X[6]*X[3];

    // 11. LacY transcription                           O -> O + M_Y            (k_s1MY)
    a[10] = c[10]*X[3];

    // 12. Leak LacY transcription                      R2O -> R2O + M_Y        (k_s0MY)
    a[11] = c[11]*X[4];

    // 13. LacY translation                             M_Y -> M_Y + Y          (k_sY)
    a[12] = c[12]*X[7];

    // 14. LacY-IPTG_ex association                     Y + Iex -> YIex         (k_p)
    a[13] = c[13]*I_ex*X[8];

    // 15. LacY-IPTG_ex dissociation                    YIex -> Y + Iex         (k_-p)
    a[14] = c[14]*X[9];

    // 16. IPTG-facilitated transport                   YIex -> Y + I           (k_ft)
    a[15] = c[15]*X[9];

    // 17. IPTG passive diffusion                       Iex -> I                (k_t)
    a[16] = V*N_A*c[16]*I_ex;

    // 18. IPTG passive diffusion                       I -> Iex                (k_t)
    a[17] = c[17]*X[5];

    // 19. LacI mRNA degradation                        M_R -> 0                (l_MR)
    a[18] = c[18]*X[0];

    // 20. LacY mRNA degradation                        MY -> 0                 (l_MY)
    a[19] = c[19]*X[7];
    
    // 21. Repressor monomer degradation                R -> 0                  (l_R)
    a[20] = c[20]*X[1];

    // 22. LacI dimer degradation                       R2 -> 0                 (l_R2)
    a[21] = c[21]*X[2];

    // 23. LacY degradation                             Y -> 0                  (l_Y)
    a[22] = c[22]*X[8];

    // 24. LacY-inducer degradation                     YIex -> I               (l_YIex)
    a[23] = c[23]*X[9];

    // 25. Repressor-inducer degradation                I2R2 -> 2I              (l_I2R2)
    a[24] = c[24]*X[6];

    // M_R, R, R2, O, R2O, I, I2R2, M_Y, Y, YIex
    // 0    1   2  3   4   5   6     7   8    9

    return a;
}
