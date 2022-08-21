#ifndef DENOISING_H
#define DENOISING_H

#include <complex>
#include <vector>
#include <math.h>

#define IS 0.25     //longitud del segmento principal sin voz
#define PI 3.141592653589793
#define wind 256


using std::vector;

class denoising
{
    public:
        /** Methods **/
        denoising();
        virtual ~denoising();
        vector<float> Weina_Im (vector<float> signal, int inc, float alpha);

        /** Variables **/
        vector<float> fsignal;  //Vector para señal filtrada
        int fs;
        int wlen;
        int inc;

    protected:

    private:
        /** Methods **/

        vector<vector<float>> enframe(vector<float> x, vector<float> win, int inc);
        vector<vector<std::complex<float>>> rfft(vector<vector<float>> y);
        vector<vector<float>> ifft(vector<vector<std::complex<float>>> Y);
        vector<vector<float>> fft_phase(vector<vector<std::complex<float>>> Y_fft);
        vector<vector<float>> fft_abs(vector<vector<std::complex<float>>> Y_fft);
        vector<float> mean(vector<vector<float>> x, int NIS);
        //std::vector<std::vector<float>> repmat;
        vector<float> hamming(int W);
        vector<float> flipframe(vector<vector<float>> x, vector<float> win, int inc);

        /** Variables **/

        int NIS;        // Número de fotogramas sin segmento inicial

};

#endif // DENOISING_H
