#include <iostream>
#include "denoising.h"
#include <complex>
#include <vector>
#include <math.h>
#include <fstream>

using std::vector;

std::ofstream file1;

//Constructor
denoising::denoising()
{

}

/** @brief (one liner)
  *
  * (documentation goes here)
  */
vector<float> denoising::flipframe(vector<vector<float>>x, vector<float>win, int inc)
{
    int nf=x[0].size();
    int len=x.size();
    int nx=(nf-1)*inc + len;
    std::cout << "nf"<<nf<<'\n';
    std::cout << "len"<<len<<'\n';

    vector<float>sig;
    vector<float>v_transition;
    vector<vector<float>>m_transition;

    for (size_t i=0; i<nf; i++){
        for (size_t j=0; j<len; j++) {
             v_transition.push_back(x[i][j]/win[j]);
        }
        m_transition.push_back(v_transition);
        v_transition.clear();
    }

    std::cout << "  -Filas    : "<< m_transition.size() <<'\n';
    std::cout << "  -Columnas : "<< m_transition[0].size()<<'\n';

    for (size_t i=0; i<nx; i++){
        sig.push_back(0);
    }

    int start=0;
    for(size_t i=0; i<nf; i++){
        for(size_t j=0; j<len; j++){
            start=(i-1)*inc;
            sig[start + j] = sig[start + j] + m_transition[i][j];
        }
    }

    std::cout << "tamano de sig: "<< sig.size() <<'\n';

    return sig;
}

/** @brief (one liner)
  *
  * (documentation goes here)
  */
vector<vector<float>> denoising::ifft(vector<vector<std::complex<float>>>Y)
{
    vector<float> x;          //para la Transformada final de un segmento
    vector<vector<float>> y;  //para la Transformada final de toda la matriz
    int N1 = Y.size();                      //cantidad de puntos del segmento
    int nfft = wind;
    std::complex<float> twiddle(0.0,0.0);   //factor Twiddle con dependencia de r
    std::complex<float> Wk(0.0,0.0);        //factor Twiddle sin dependencia de r
    std::complex<float> fe(0.0,0.0);        //Transformada de muestras pares de un segmento
    std::complex<float> fo(0.0,0.0);        //Transformada de muestras impares de un segmento

    for(size_t i=0; i<Y.size();i++){
        for(size_t n=0; n< nfft; n++){
            for(size_t r=0; r <= N1/2-1; r++){
                twiddle=(std::complex<float>(cos(-2*PI*n*r/(N1/2)),sin(-2*PI*n*r/(N1/2))));
                fe+= Y[2*r][i]*twiddle;
                fo+= Y[2*r+1][i]*twiddle;
            }
            Wk=std::complex<float> (cos(-2*PI*n/N1),sin(-2*PI*n/N1));
            x.push_back((1/(float)N1)*(std::real(fe+Wk*fo)));
            fe=std::complex<float>(0.0,0.0);
            fo=std::complex<float>(0.0,0.0);

        }
        y.push_back(x);
        x.clear();
    }

    std::cout << "  -Y[0]    : "<< Y[0].size() <<'\n';
    std::cout << "  -Filas    : "<< y.size() <<'\n';
    std::cout << "  -Columnas : "<< y[0].size()<<'\n';


    return y;
}

/** @brief (one liner)
  *
  * (documentation goes here)
  */
vector <float> denoising::mean(vector<vector<float>> x, int NIS)
{
    float N_s{0.0};
    vector<float> Ns;

    for (size_t i=0; i<x[0].size(); i++){
        for (size_t j=0; j<NIS; j++) {
                N_s += x[j][i];
        }
        Ns.push_back(N_s/(NIS)); //promedio
        N_s=0.0;
    }

    file1.open("C:/CodeB/test1/noise.txt");
    for (size_t j=0; j < Ns.size(); j++){
        file1 << Ns[j]<<'\n';
    }
    file1.close();

    return Ns;


}

/** @brief (one liner)
  *
  * (documentation goes here)
  */
vector<vector<float>> denoising::fft_abs(vector<vector<std::complex<float>>> Y_fft)
{
    vector<float>ROWf;
    vector<vector<float>>Yabs;
    for (size_t i=0; i<Y_fft.size(); i++){
        for (size_t j=0; j < Y_fft[i].size(); j++) {
            ROWf.push_back(std::abs(Y_fft[i][j]));
        }
        Yabs.push_back(ROWf);
        ROWf.clear();
    }

    return Yabs;
}

/** @brief (one liner)
  *
  * (documentation goes here)
  */
vector<vector<float>> denoising::fft_phase(vector<vector<std::complex<float>>> Y_fft)
{
    std::vector<float>ROWf;
    vector<vector<float>> Yphase;
    for (size_t i=0; i<Y_fft.size(); i++){
        for (size_t j=0; j < Y_fft[i].size(); j++) {
            ROWf.push_back(std::arg(Y_fft[i][j]));
        }
        Yphase.push_back(ROWf);

        ROWf.clear();
    }

    return Yphase;
}

/** @brief (one liner)
  *
  * (documentation goes here)
  */
vector<vector<std::complex<float>>> denoising::rfft(vector<vector<float>> y)
{
    vector<std::complex<float>> X;          //para la Transformada final de un segmento
    vector<vector<std::complex<float>>> Y;  //para la Transformada final de toda la matriz
    int N1 = y.size();                      //cantidad de puntos del segmento
    int nfft = wind;
    std::complex<float> twiddle(0.0,0.0);   //factor Twiddle con dependencia de r
    std::complex<float> Wn(0.0,0.0);        //factor Twiddle sin dependencia de r
    std::complex<float> fe(0.0,0.0);        //Transformada de muestras pares de un segmento
    std::complex<float> fo(0.0,0.0);        //Transformada de muestras impares de un segmento

    for(size_t i=0; i<y[0].size();i++){
        for(size_t k=0; k< nfft; k++){
          if (k< nfft/2+1){ //primera mitad
            for(size_t r=0; r <= N1/2-1; r++){
                twiddle=(std::complex<float>(cos(-2*PI*k*r/(N1/2)),sin(-2*PI*k*r/(N1/2))));
                fe+= y[2*r][i]*twiddle;
                fo+= y[2*r+1][i]*twiddle;
            }
            Wn=std::complex<float> (cos(-2*PI*k/N1),sin(-2*PI*k/N1));
            X.push_back(fe + Wn*fo);
            fe=std::complex<float>(0.0,0.0);
            fo=std::complex<float>(0.0,0.0);
          }
          else{//conjugado
              X.push_back(std::conj(X[nfft-k]));
          }
        }
        Y.push_back(X);
        X.clear();
    }

    return Y;

}

/** @brief (one liner)
  *
  * (documentation goes here)
  */
vector<vector<float>> denoising::enframe(vector<float> x, vector<float> win, int inc)
{
    int nx=x.size();
    int nwin=win.size();
    int nli=nx-nwin+inc;
    int nf=floor(nli/inc);

    vector<vector<float>>f;
    vector<float> v_transition;

    for (size_t i=0; i<nwin; i++){
      for (size_t j=0; j<nf; j++){
        v_transition.push_back(x[inc*j+i]*win[i]);
      }
      f.push_back(v_transition);
      v_transition.clear();
    }



    return f;

}

/** @brief (one liner)
  *
  * (documentation goes here)
  */
vector<float> denoising::hamming(int W)
{
    vector<float> win_hamming;
    for(size_t i=0; i < W; i++){
        win_hamming.push_back(0.54 - 0.46*cos(2*PI*(float(i)/float(W-1))));
    }

    return win_hamming;
}


/**
  *
  * (documentation goes here)
  */
vector<float> denoising::Weina_Im(vector<float> signal, int inc, float alpha)
{
    int Length;
    int framesize;
    int framenum;
    vector<float> wnd;
    vector<vector<float>> y;
    vector<vector<std::complex<float>>> Y_fft;
    vector<vector<float>> Y_phase;
    vector<vector<float>> Y_abs;
    vector<vector<float>> Y_fft2;
    vector<float> v_transition;
    vector<float> noise;
    vector<float>snr_x;
    vector<float>snr_x_q;
    vector<float>Hw;
    vector<float>M;
    vector<std::complex<float>>Mn;
    vector<vector<std::complex<float>>>SigF;
    vector<vector<float>>sigT;

    NIS = floor((IS*fs-wlen)/inc+1);

    Length = signal.size();

    framesize = wind;
    wnd = hamming(framesize);

    std::cout<< "Enmarcando Senal... ";
    y = enframe(signal,wnd,inc);
    framenum = y[0].size();
    std::cout<< "Finalizado! - Numero de marcos: "<<framenum<<'\n';

    std::cout<< "Calculando FFT... ";
    Y_fft = rfft(y);
    std::cout<< "Finalizado!"<<'\n';

    std::cout<< "Calculando Magnitud de FFT... ";
    Y_abs = fft_abs(Y_fft);
    std::cout<< "Finalizado!"<<'\n';

    std::cout<< "Calculando Fase de FFT... ";
    Y_phase = fft_phase(Y_fft);
    std::cout<< "Finalizado!"<<'\n';

    std::cout<< "Calculando potencia... ";
    for (size_t i=0; i<Y_abs.size(); i++){
        for (size_t j=0; j < Y_abs[i].size(); j++) {
            v_transition.push_back(pow(Y_abs[i][j],2));
        }
        Y_fft2.push_back(v_transition);
        v_transition.clear();
    }
    std::cout<< "Finalizado!"<<'\n';

    std::cout<< "Calculando Ruido... ";
    noise=mean(Y_fft2, NIS);
    std::cout<< "Finalizado!"<<'\n';

    std::cout<< "Algoritmo de Weiner... ";
    v_transition.clear();
    snr_x_q.push_back(0.96);
    for (size_t i=0; i < framenum; i++){

        if(i==0){
            for (size_t j=0; j < framesize; j++){

            if((Y_fft2[i][j]/noise[j])-1.0>0.0){
                snr_x.push_back(alpha*0.95+ (1.0-alpha)*((Y_fft2[i][j]/noise[j])-1.0));
            }else{
                snr_x.push_back(alpha*0.95);
            }
            Hw.push_back(snr_x[j]/(1.0+snr_x[j]));
            M.push_back(Y_abs[i][j]*Hw[j]);
            Mn.push_back(std::complex<float>(M[j]*cos(Y_phase[i][j]),M[j]*sin(Y_phase[i][j])));
            snr_x_q.push_back((pow(M[j],2))/noise[j]);
            }

            SigF.push_back(Mn);

            snr_x.clear();
            Hw.clear();
            M.clear();
            Mn.clear();
        }
        else{
            for (size_t j=0; j < framesize; j++){

                if((Y_fft2[i][j]/noise[j])-1.0>0.0){
                    snr_x.push_back(alpha*snr_x_q[j]+ (1.0-alpha)*((Y_fft2[i][j]/noise[j])-1.0));
                }else{
                    snr_x.push_back(alpha*snr_x_q[j]);
                }

                Hw.push_back(snr_x[j]/(1.0+snr_x[j]));
                M.push_back(Y_abs[i][j]*Hw[j]);
                Mn.push_back(std::complex<float>(M[j]*cos(Y_phase[i][j]),M[j]*sin(Y_phase[i][j])));
                snr_x_q.clear();
                snr_x_q.push_back((pow(M[j],2))/noise[j]);
            }
            SigF.push_back(Mn);

            snr_x.clear();
            Hw.clear();
            M.clear();
            Mn.clear();

        }

    }
    /*
        file1.open("C:/CodeB/test1/sig.txt");
        for (size_t j=0; j < SigF[0].size(); j++){
            file1 << SigF[0][j]<<'\n';
        }
        file1.close();

        file1.open("C:/CodeB/test1/snr_x_q.txt");
        for (size_t j=0; j < snr_x_q.size(); j++){
            file1 << snr_x_q[j]<<'\n';
        }
        file1.close();
        std::cout << "SigF. Tamano:"<<'\n';
        std::cout << "  -Filas    : "<< SigF.size() <<'\n';
        std::cout << "  -Columnas : "<< SigF[0].size()<<'\n';
*/

        std::cout<< "Finalizado!"<<'\n';

        std::cout << "Volviendo al dominio del tiempo... ";
        sigT = ifft(SigF);
        std::cout<< "Finalizado!"<<'\n';

        std::cout << "Reconstruyendo senal... ";
        //fsignal = flipframe(sigT,wnd,inc);
        std::cout<< "Finalizado!"<<'\n';
    return fsignal;
}



//Destructor
denoising::~denoising()
{
    //dtor
}
