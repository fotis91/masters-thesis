/*
 * Author: Fotios Kostarelos
 *
 * This HLS IP in C code, calculates the AmplitudeAsymmetry  
 * between two EEGG channels as described in the thesis.
 * 
 *
 *   
 */
#include <ap_axi_sdata.h>
#include <hls_stream.h>
#include <math.h>
#define FS 125 // sampling frequency
#define M 7 /* Number of FFT Stages = Log2FFTpoints, FFTpoints=128, log2(128)=7 */
#define N 8192 //input length
#define SIZE 128 //FFT SIZE
#define SIZE2 SIZE>>1	/* SIZE/2 */
#define Navg 127 // number of iterations for spectrum calculation as defined by Welch's method(2*N/SIZE-1)
#define samplesread SIZE-1

#define betaBandStart 11
#define betaBandEnd 21

#define AlphaBandStart 6
#define AlphaBandEnd 12

//This structure defines what signals are implemented for the custom streaming variable "varvara" from the standard AXI4-Stream signals
typedef struct{
   float data;
   ap_uint<1> last;
}varvara;

float hamming_window[SIZE] = {0.1273947836455926,0.1282910832476171,0.1309777886583703,0.1354483250594968,0.1416917522992888,0.1496927916650700,0.1594318632727663,0.1708851339821634,0.1840245757205936,0.1988180340723275,0.2152293069658165,0.2332182332662292,0.2527407910564775,0.2737492053662246,0.2961920650852438,0.3200144487750193,0.3451580590707104,0.3715613653445714,0.3991597542817092,0.4278856879996919,0.4576688693250643,0.4884364138223117,0.5201130281542863,0.5526211943376206,0.5858813594422227,0.6198121302706268,0.6543304725407861,0.6893519140848757,0.7247907515668432,0.7605602602128315,0.7965729060412342,0.8327405600730173,0.8689747139980977,0.9051866967700101,0.9412878915988160,0.9771899528112388,1.0128050220473304,1.0480459432645999,1.0828264760234558,1.1170615065320122,1.1506672559337996,1.1835614853286665,1.2156636970251495,1.2468953315318121,1.2771799598054854,1.3064434702859438,1.3346142502593119,1.3616233611063759,1.3874047070069355,1.4118951966873552,1.4350348978154785,1.4567671836650915,1.4770388716910086,1.4958003536756734,1.5130057171287807,1.5286128576428328,1.5425835819296794,1.5548837012858960,1.5654831152582693,1.5743558853046491,1.5814802982699081,1.5868389195216730,1.5904186356157881,1.5922106863871124,1.5922106863871124,1.5904186356157881,1.5868389195216730,1.5814802982699083,1.5743558853046491,1.5654831152582693,1.5548837012858965,1.5425835819296796,1.5286128576428326,1.5130057171287810,1.4958003536756734,1.4770388716910086,1.4567671836650915,1.4350348978154785,1.4118951966873554,1.3874047070069360,1.3616233611063757,1.3346142502593121,1.3064434702859440,1.2771799598054856,1.2468953315318121,1.2156636970251495,1.1835614853286669,1.1506672559337998,1.1170615065320124,1.0828264760234563,1.0480459432646003,1.0128050220473304,0.9771899528112393,0.9412878915988165,0.9051866967700100,0.8689747139980981,0.8327405600730177,0.7965729060412345,0.7605602602128316,0.7247907515668432,0.6893519140848757,0.6543304725407858,0.6198121302706271,0.5858813594422231,0.5526211943376202,0.5201130281542864,0.4884364138223117,0.4576688693250643,0.4278856879996923,0.3991597542817096,0.3715613653445712,0.3451580590707106,0.3200144487750195,0.2961920650852439,0.2737492053662250,0.2527407910564773,0.2332182332662291,0.2152293069658167,0.1988180340723276,0.1840245757205937,0.1708851339821637,0.1594318632727663,0.1496927916650699,0.1416917522992889,0.1354483250594969,0.1309777886583704,0.1282910832476171,0.1273947836455926};
float sinearray[SIZE2]={0.000000000000000000,  -0.049067676067352295,  -0.098017141222953796,  -0.146730467677116390,  -0.195090323686599730,  -0.242980197072029110,  -0.290284663438797000,  -0.336889833211898800,  -0.382683426141738890,  -0.427555054426193240,  -0.471396714448928830,  -0.514102756977081300,  -0.555570244789123540,  -0.595699310302734380,  -0.634393334388732910,  -0.671558976173400880,  -0.707106828689575200,  -0.740951180458068850,  -0.773010551929473880,  -0.803207635879516600,  -0.831469714641571040,  -0.857728719711303710,  -0.881921350955963130,  -0.903989374637603760,  -0.923879623413085940,  -0.941544175148010250,  -0.956940412521362300,  -0.970031321048736570,  -0.980785369873046880,  -0.989176571369171140,  -0.995184779167175290,  -0.998795449733734130,  -1.000000000000000000,  -0.998795449733734130,  -0.995184659957885740,  -0.989176452159881590,  -0.980785191059112550,  -0.970031142234802250,  -0.956940174102783200,  -0.941543877124786380,  -0.923879325389862060,  -0.903989076614379880,  -0.881921112537384030,  -0.857728481292724610,  -0.831469535827636720,  -0.803207516670227050,  -0.773010492324829100,  -0.740951240062713620,  -0.707106947898864750,  -0.671559214591979980,  -0.634393632411956790,  -0.595699727535247800,  -0.555570781230926510,  -0.514103353023529050,  -0.471397459506988530,  -0.427555918693542480,  -0.382684379816055300,  -0.336890906095504760,  -0.290285855531692500,  -0.242981463670730590,  -0.195091724395751950,  -0.146731987595558170,  -0.098018750548362732,  -0.049069393426179886 };
float cosinearray[SIZE2]={1.000000000000000000,  0.998795449733734130,  0.995184719562530520,  0.989176511764526370,  0.980785250663757320,  0.970031261444091800,  0.956940352916717530,  0.941544055938720700,  0.923879563808441160,  0.903989315032958980,  0.881921291351318360,  0.857728600502014160,  0.831469595432281490,  0.803207516670227050,  0.773010432720184330,  0.740951061248779300,  0.707106709480285640,  0.671558856964111330,  0.634393215179443360,  0.595699191093444820,  0.555570125579833980,  0.514102578163146970,  0.471396565437316890,  0.427554905414581300,  0.382683217525482180,  0.336889594793319700,  0.290284395217895510,  0.242979884147644040,  0.195090010762214660,  0.146730139851570130,  0.098016783595085144,  0.049067292362451553,  -0.000000401339264045,  -0.049068093299865723,  -0.098017580807209015,  -0.146730929613113400,  -0.195090800523757930,  -0.242980659008026120,  -0.290285170078277590,  -0.336890369653701780,  -0.382683962583541870,  -0.427555501461029050,  -0.471397042274475100,  -0.514102995395660400,  -0.555570363998413090,  -0.595699369907379150,  -0.634393274784088130,  -0.671558856964111330,  -0.707106590270996090,  -0.740950882434844970,  -0.773010194301605220,  -0.803207218647003170,  -0.831469237804412840,  -0.857728242874145510,  -0.881920874118804930,  -0.903988897800445560,  -0.923879146575927730,  -0.941543698310852050,  -0.956939995288848880,  -0.970030903816223140,  -0.980785012245178220,  -0.989176273345947270,  -0.995184540748596190,  -0.998795390129089360 };

unsigned int reverse_bits(unsigned int input);
void bit_reverse(float X_R[SIZE], float X_I[SIZE]);
void fft(float X_R[SIZE], float X_I[SIZE]);

void AmplitudeAsymmetry(hls::stream<varvara>& a,hls::stream<varvara>& b,float * AA1, float * AA2){

	//variables declaration
	int i=0,j=0,k=0,l=0,m=0,n=0;
	float x1[SIZE],y1[SIZE],x1i[SIZE],y1i[SIZE],areadstream[SIZE/2],breadstream[SIZE/2];
	float magsmpax[SIZE/2],magsmpay[SIZE/2],Phzx[SIZE/2],Phzy[SIZE/2];
    float AAbetaalpha[2],betaAAacc=0,alphaAAacc=0,alphaAA=0,betaAA=0;
    float suma=0,sumb=0;
    varvara tmpa,tmpb;
    varvara tmpAA1,tmpAA2;

    for(n=0;n<SIZE/2;n++){
    	Phzx[n]=0;//set arrays values to zero to avoid reading trash
    	Phzy[n]=0;//set arrays values to zero to avoid reading trash
    }

	for(i=0;i<Navg;i++){//number of iterations for periodograms(welch's method)


		for(j=0;j<=samplesread;j++){ //read 128 sample
					if(i==0){//read the fist 128 samples
						tmpa=a.read();
						tmpb=b.read();
						if(j>samplesread>>1){//store in the buffer the second half of them for use in the next iteration
							areadstream[j-(samplesread>>1)-1] = tmpa.data;
							breadstream[j-(samplesread>>1)-1] = tmpb.data;
						}
					}
					else{
						if(j<=(samplesread>>1)){//read the previously stored values from the buffer
							tmpa.data = areadstream[j];
							tmpb.data = breadstream[j];
						}
						else{//store in the buffer the second half of them for use in the next iteration
							tmpa=a.read();
							tmpb=b.read();
							areadstream[j-(samplesread>>1)-1] = tmpa.data;
							breadstream[j-(samplesread>>1)-1] = tmpb.data;
						}

					}

			/////////// Windowing
			x1[j] = tmpa.data*hamming_window[j];//read samples and multiply with the window function
			y1[j] = tmpb.data*hamming_window[j];//read samples and multiply with the window function
			x1i[j] = 0;//set imaginary part equal to zero for x
			y1i[j] = 0;//set imaginary part equal to zero for x


		}
			/////////// FFT
			fft(x1,x1i);
			fft(y1,y1i);

			// spectral calculations
			 for(k=0;k<SIZE/2;k++){
				 magsmpax[k] = (x1[k]*x1[k]+x1i[k]*x1i[k])/Navg;//squaring and averaging
				 magsmpay[k] = (y1[k]*y1[k]+y1i[k]*y1i[k])/Navg;//squaring and averaging
				 if(k==0){
				  	Phzx[k] = Phzx[k] + (magsmpax[k])/(SIZE*FS);//normalize the values
				   	Phzy[k] = Phzy[k] + (magsmpay[k])/(SIZE*FS);//normalize the values
				  	}
				 else{
				   	Phzx[k] = Phzx[k] + (2*magsmpax[k])/(SIZE*FS);//normalize and fold the values except of the DC component
				   	Phzy[k] = Phzy[k] + (2*magsmpay[k])/(SIZE*FS);//normalize and fold the values except of the DC component
				   	}

		        }


	}

		// feature calculations
		for(l=AlphaBandStart;l<=betaBandEnd;l++){

			if(l<=AlphaBandEnd){
				alphaAA = (sqrt(Phzx[l])-sqrt(Phzy[l]))/(sqrt(Phzx[l])+sqrt(Phzy[l]));
				suma = suma + alphaAA;
			}

			if(l>=betaBandStart){
               	betaAA = (sqrt(Phzx[l])-sqrt(Phzy[l]))/(sqrt(Phzx[l])+sqrt(Phzy[l]));
	        	sumb = sumb + betaAA;}
	        }

			alphaAA=suma/(AlphaBandEnd-AlphaBandStart+1);
			betaAA=sumb/(betaBandEnd-betaBandStart+1);


			//set the output, when the last signal goes high
		    if( tmpb.last)
		    {
		    *AA1 = alphaAA;
			*AA2 = betaAA;
		    }

}

unsigned int reverse_bits(unsigned int input) {
	int i, rev = 0;
	for (i = 0; i < M; i++) {
		rev = (rev << 1) | (input & 1);
		input = input >> 1;
	}
	return rev;
}

void bit_reverse(float X_R[SIZE], float X_I[SIZE]) {
	unsigned int reversed;
	unsigned int i;
	float temp;

	for (i = 0; i < SIZE; i++) {
		reversed = reverse_bits(i); // Find the bit reversed index
		if (i <= reversed) {
			// Swap the real values
			temp = X_R[i];
			X_R[i] = X_R[reversed];
			X_R[reversed] = temp;

			// Swap the imaginary values
			temp = X_I[i];
			X_I[i] = X_I[reversed];
			X_I[reversed] = temp;
		}
	}
}



void fft(float X_R[SIZE], float X_I[SIZE]) {
	float temp_R; // temporary storage complex variable
	float temp_I; // temporary storage complex variable
	int i, j, k;	// loop indexes
	int i_lower;	// Index of lower point in butterfly
	int step, stage, DFTpts,initposition,position;
	int numBF;			// Butterfly Width
	int N2 = SIZE2; // N2=N>>1

	bit_reverse(X_R, X_I);

	step = N2;
	float a, e, c, s, pioverangle ,pos;

stage_loop:
	for (stage = 1; stage <= M; stage++) { // Do M stages of butterflies
		DFTpts = 1 << stage;								 // DFT = 2^stage = points in sub DFT
		numBF = DFTpts >> 1;//numBF = DFTpts / 2;									 // Butterfly WIDTHS in sub-DFT
		//k = 0;
		//e = unitcirclesegmentsarray[stage-1];//e = -6.283185307178 / DFTpts;
		//a = 0.0;
		initposition=(SIZE/2)>>stage-1;
		position=0;
	// Perform butterflies for j-th stage
	butterfly_loop:
		for (j = 0; j < numBF; j++) {

		    if(position==0){
		    c=cosinearray[0];
		    s=sinearray[0];
		    		}
		    else{
		    c = cosinearray[position];
		    s = sinearray[position];
		    		}

		    position = position + initposition;

			//a = a + e;
		// Compute butterflies that use same W**k
		dft_loop:
			for (i = j; i < SIZE; i += DFTpts) {
				i_lower = i + numBF; // index of lower point in butterfly
				temp_R = X_R[i_lower] * c - X_I[i_lower] * s;
				temp_I = X_I[i_lower] * c + X_R[i_lower] * s;
				X_R[i_lower] = X_R[i] - temp_R;
				X_I[i_lower] = X_I[i] - temp_I;
				X_R[i] = X_R[i] + temp_R;
				X_I[i] = X_I[i] + temp_I;
			}
			//k += step;
		}
		//step = step / 2;
	}
}
