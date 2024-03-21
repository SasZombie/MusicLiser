#include <iostream>
#include <raylib.h>
#include <vector>
#include <climits>
#include <cstring>
#include <complex>

using namespace std::complex_literals;

struct Frame{
    float left, right;
};

void callback(void *buffedData, unsigned int frames);
void fft(double in2[], size_t stride, std::complex<double>out2[], size_t n);
float amp(std::complex<double> z);


std::vector<Frame> globalFrames;
unsigned int globalFramesCount = 1;
//const int N = 256;
const int N = (1<<10);

float max_amp;
double in[N];
std::complex<double> out[N];



void fft(double in2[], size_t stride, std::complex<double>out2[], size_t n)
{
    if( n == 1)
    {
        out2[0] = in2[0];
        return;
    }
    fft(in2, stride*2, out2, n/2);

    fft(in2 + stride, stride*2, out2 + n/2, n/2);

    for (size_t k = 0; k < n/2; ++k)
    {
        double t = (double) k/n;
        std::complex<double> v = std::exp((-2) * PI * t * 1i) * out2[k+n/2];
        std::complex<double> e = out2[k];
        out2[k] = e + v;
        out2[k+n/2] = e-v;
    } 

}


float amp(std::complex<double> z)
{
        float c = std::abs(z.real());
        float d = std::abs(z.imag());

        if(c < d)
            return d;
        return c;
}


void callback(void *bufferData, unsigned int frames)
{

    globalFramesCount = frames;
    globalFrames.resize(frames);
    std::memcpy(globalFrames.data(), bufferData, frames * sizeof(Frame));


    
    Frame *fs = (Frame*)bufferData;

    for (size_t i = 0; i < frames; ++i)
    {
        memmove(in, in+1, (N-1)*sizeof(in[0]));

        in[i] = fs[i].left;
    }

    fft(in, 1, out, N);

    max_amp = 0;

    for (size_t i = 0; i < frames; ++i)
    {
        float a = amp(out[i]);
        if(max_amp < a) max_amp = a; 
    }
    
   
  
}


int main()
{
    InitWindow(800, 600, "Musializer");
    SetTargetFPS(60);
    InitAudioDevice();
    Music music = LoadMusicStream("/home/saszombie/Downloads/Into You (Sped Up).mp3");

    PlayMusicStream(music);
    AttachAudioStreamProcessor(music.stream, callback);
    SetMusicVolume(music, 0.1f);
    int width = GetRenderWidth();
    int height = GetRenderHeight();

    while (!WindowShouldClose())
    {

        if(IsFileDropped())
        {
            FilePathList draoppedFiles = LoadDroppedFiles();
            
            if(draoppedFiles.count > 0)
            {
                music = LoadMusicStream(draoppedFiles.paths[0]);
                PlayMusicStream(music);
                SetMusicVolume(music, 0.1f);
                AttachAudioStreamProcessor(music.stream, callback);

            }

            UnloadDroppedFiles(draoppedFiles);
        }
        UpdateMusicStream(music);

        BeginDrawing();
        ClearBackground(CLITERAL(Color){0x18, 0x18, 0x18, 0xFF});
        
         
        float step = 1.06f;
        size_t m = 0;

        for (float i = 20.0f; static_cast<size_t>(i) < N; i*=step)
        {
            ++m;
        }
        
        
        float cell_width = static_cast<float>(width)/m;
        m = 0;

        for(float f = 20.0f; static_cast<size_t>(f) < N; f*=step)
        {
            float f1 = f*step;
            float a = 0.0f;

            for (size_t q = static_cast<size_t>(f); q < N && q < static_cast<size_t>(f1); ++q)
                a+=amp(out[q]);

            a/= static_cast<size_t>(f1) - static_cast<size_t>(f+1);

            float t = a/max_amp;
            
            DrawRectangle(m*cell_width, height/2 - height/2 * t, cell_width , height/2*t, PURPLE); 
            //DrawLine(m*cell_width, height/2-height/2 * t, cell_width, height/2*t, PURPLE);
            //DrawCircle(m*cell_width, height/2-height/2 * t, 2   , PURPLE);
        
            m+=1;

            for(unsigned long q = 1; q < globalFramesCount; ++q)
            {
                if(globalFrames.size() <= globalFramesCount)
                {
                    if(globalFrames.size() > q)
                    { 
                        float t2 = globalFrames[q].left;
                        if(t2 > 0)
                        { 
                            DrawRectangle(q*cell_width, (height- height/2 * t2),1 , height/2*t2, RED);
                        }
                        else
                        {
                            DrawRectangle(q*cell_width, height, 1, height/2*t2, RED);
                        }
                    }
                }
            }
        }
   


        EndDrawing();
    }
    
}